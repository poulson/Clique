/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "clique.hpp"

// This is the part of the symbolic factorization that requires fine-grain 
// parallelism: we assume that the upper floor(log2(commSize)) levels of the
// tree are balanced.
//
// TODO: Generalize so that the depth can be less than or equal to 
// floor(log2(commSize)). This would allow for the use of more processes in the
// factorization.
//
// TODO: Generalize to support more than just a power-of-two number of 
//       processes. This should be relatively straightforward.
void clique::symbolic::DistSymmetricFactorization
( const DistSymmOrig&  distOrig,
  const LocalSymmFact& localFact, DistSymmFact&  distFact, 
        bool storeFactRecvIndices )
{
#ifndef RELEASE
    PushCallStack("symbolic::DistSymmetricFactorization");
#endif
    const unsigned numSupernodes = distOrig.supernodes.size();
    distFact.supernodes.resize( numSupernodes );
    if( numSupernodes == 0 )
        return;

    const unsigned commRank = mpi::CommRank( distOrig.comm );
    const unsigned commSize = mpi::CommSize( distOrig.comm );
#ifndef RELEASE
    // Use the naive algorithm for finding floor(log2(commSize)) since it
    // will be an ignorable portion of the overhead and will be more easily
    // extended to more general integer types.
    unsigned temp = commSize;
    unsigned log2CommSize = 0;
    while( temp >>= 1 )
        ++log2CommSize;
    if( log2CommSize != numSupernodes )
        throw std::runtime_error("Invalid distributed tree depth");
    if( 1u<<log2CommSize != commSize )
        throw std::runtime_error
        ("Power-of-two number of procs currently required");
#endif

    // The bottom node is already computed, so just copy it over
    const LocalSymmFactSupernode& topLocalSN = localFact.supernodes.back();
    DistSymmFactSupernode& bottomDistSN = distFact.supernodes[0];
    mpi::CommSplit( distOrig.comm, commRank, 0, bottomDistSN.comm );
    unsigned gridHeight = 
        static_cast<unsigned>(sqrt(static_cast<double>(commSize)));
    while( commSize % gridHeight != 0 )
        ++gridHeight;
    bottomDistSN.gridHeight = gridHeight;
    bottomDistSN.size = topLocalSN.size;
    bottomDistSN.offset = topLocalSN.offset;
    bottomDistSN.lowerStruct = topLocalSN.lowerStruct;
    bottomDistSN.origLowerRelIndices = topLocalSN.origLowerRelIndices;
    bottomDistSN.leftChildRelIndices = topLocalSN.leftChildRelIndices;
    bottomDistSN.rightChildRelIndices = topLocalSN.rightChildRelIndices;
    bottomDistSN.leftChildSize = -1; // not needed, could compute though
    bottomDistSN.rightChildSize = -1; // not needed, could compute though
    bottomDistSN.leftChildFactColIndices.clear();
    bottomDistSN.leftChildFactRowIndices.clear();
    bottomDistSN.rightChildFactColIndices.clear();
    bottomDistSN.rightChildFactRowIndices.clear();
    bottomDistSN.numChildFactSendIndices.clear();
    bottomDistSN.childFactRecvIndices.clear();

    // Perform the distributed part of the symbolic factorization
    std::vector<int>::iterator it;
    std::vector<int> sendBuffer, recvBuffer;
    std::vector<int> childrenStruct, partialStruct, fullStruct,
                    supernodeIndices;
    for( unsigned k=1; k<numSupernodes; ++k )
    {
        const DistSymmOrigSupernode& origSN = distOrig.supernodes[k];
        const DistSymmFactSupernode& factChildSN = distFact.supernodes[k-1];
        DistSymmFactSupernode& factSN = distFact.supernodes[k];
        factSN.size = origSN.size;
        factSN.offset = origSN.offset;

        // Determine our partner based upon the bits of 'commRank'
        const unsigned powerOfTwo = 1u << (k-1);
        const unsigned partner = commRank ^ powerOfTwo; // flip the k-1'th bit
        const bool onLeft = (commRank & powerOfTwo) == 0; // check k-1'th bit 

        // Create this level's communicator
        const unsigned teamSize  = powerOfTwo;
        const int teamColor = commRank & !(powerOfTwo-1);
        const int teamRank  = commRank &  (powerOfTwo-1);
        mpi::CommSplit( distOrig.comm, teamColor, teamRank, factSN.comm );
        gridHeight = static_cast<unsigned>(sqrt(static_cast<double>(teamSize)));
        while( teamSize % gridHeight != 0 )
            ++gridHeight;
        factSN.gridHeight = gridHeight;
        const unsigned gridWidth = teamSize / gridHeight;
        const unsigned gridRow = teamRank % gridHeight;
        const unsigned gridCol = teamRank / gridHeight;

        // Retrieve the child grid information
        const unsigned childTeamRank = mpi::CommRank( factChildSN.comm );
        const unsigned childTeamSize = mpi::CommSize( factChildSN.comm );
        const unsigned childGridHeight = factChildSN.gridHeight;
        const unsigned childGridWidth = childTeamSize / childGridHeight;
        const unsigned childGridRow = childTeamRank % childGridHeight;
        const unsigned childGridCol = childTeamRank / childGridHeight;

        // SendRecv the message lengths
        const int myChildSize = factChildSN.size;
        const int myChildLowerStructSize = factChildSN.lowerStruct.size();
        const int initialSends[2] = { myChildSize, myChildLowerStructSize };
        int initialRecvs[2];
        mpi::SendRecv
        ( initialSends, 2, partner, 0,
          initialRecvs, 2, partner, 0, distOrig.comm );
        const int theirChildSize = initialRecvs[0];
        const int theirChildLowerStructSize = initialRecvs[1];
        // Perform the exchange
        sendBuffer.resize( myChildLowerStructSize );
        recvBuffer.resize( theirChildLowerStructSize );
        std::memcpy
        ( &sendBuffer[0], &factChildSN.lowerStruct[0], 
          myChildLowerStructSize*sizeof(int) );
        mpi::SendRecv
        ( &sendBuffer[0], myChildLowerStructSize, partner, 0,
          &recvBuffer[0], theirChildLowerStructSize, partner, 0, 
          distOrig.comm );
        
        // Union the two child lower structures
        childrenStruct.resize
        ( myChildLowerStructSize+theirChildLowerStructSize );
        it = std::set_union
        ( sendBuffer.begin(), sendBuffer.end(),
          recvBuffer.begin(), recvBuffer.end(), childrenStruct.begin() );
        const int childrenStructSize = int(it-childrenStruct.begin());
        childrenStruct.resize( childrenStructSize );

        // Union the lower structure of this supernode
        partialStruct.resize( childrenStructSize + origSN.lowerStruct.size() );
        it = std::set_union
        ( childrenStruct.begin(), childrenStruct.end(),
          origSN.lowerStruct.begin(), origSN.lowerStruct.end(),
          partialStruct.begin() );
        const int partialStructSize = int(it-partialStruct.begin());
        partialStruct.resize( partialStructSize );

        // Union again with the supernode indices
        supernodeIndices.resize( origSN.size );
        for( int i=0; i<origSN.size; ++i )
            supernodeIndices[i] = origSN.offset + i;
        fullStruct.resize( origSN.size + partialStructSize );
        it = std::set_union
        ( supernodeIndices.begin(), supernodeIndices.end(),
          partialStruct.begin(), partialStruct.end(), 
          fullStruct.begin() );
        const int fullStructSize = int(it-fullStruct.begin());
        fullStruct.resize( fullStructSize );

        // Construct the relative indices of the original lower structure
        const int numOrigLowerIndices = origSN.lowerStruct.size();
        it = fullStruct.begin();
        for( int i=0; i<numOrigLowerIndices; ++i )
        {
            const int index = origSN.lowerStruct[i];
            it = std::lower_bound( it, fullStruct.end(), index );
            factSN.origLowerRelIndices[index] = *it;
        }

        // Construct the relative indices of the children
        int numLeftIndices, numRightIndices;
        const int *leftIndices, *rightIndices;
        if( onLeft )
        {
            factSN.leftChildSize = myChildSize;
            factSN.rightChildSize = theirChildSize;
            leftIndices = &sendBuffer[0];
            rightIndices = &recvBuffer[0];
            numLeftIndices = sendBuffer.size();
            numRightIndices = recvBuffer.size();
        }
        else
        {
            factSN.leftChildSize = theirChildSize;
            factSN.rightChildSize = myChildSize;
            leftIndices = &recvBuffer[0];
            rightIndices = &sendBuffer[0];
            numLeftIndices = recvBuffer.size();
            numRightIndices = sendBuffer.size();
        }
        factSN.leftChildRelIndices.resize( numLeftIndices );
        it = fullStruct.begin();
        for( int i=0; i<numLeftIndices; ++i )
        {
            const int index = leftIndices[i];
            it = std::lower_bound( it, fullStruct.end(), index );
            factSN.leftChildRelIndices[i] = *it;
        }
        factSN.rightChildRelIndices.resize( numRightIndices );
        it = fullStruct.begin();
        for( int i=0; i<numRightIndices; ++i )
        {
            const int index = rightIndices[i];
            it = std::lower_bound( it, fullStruct.end(), index );
            factSN.rightChildRelIndices[i] = *it;
        }

        // Form lower structure of this node by removing the supernode indices
        factSN.lowerStruct.resize( fullStructSize );
        it = std::set_difference
        ( fullStruct.begin(), fullStruct.end(),
          supernodeIndices.begin(), supernodeIndices.end(),
          factSN.lowerStruct.begin() );
        factSN.lowerStruct.resize( int(it-factSN.lowerStruct.begin()) );

        // Fill numChildFactSendIndices so that we can reuse it for many facts.
        factSN.numChildFactSendIndices.resize( teamSize );
        std::memset
        ( &factSN.numChildFactSendIndices[0], 0, teamSize*sizeof(int) );
        const std::vector<int>& myChildRelIndices = 
            ( onLeft ? factSN.leftChildRelIndices 
                     : factSN.rightChildRelIndices );
        const int updateSize = factSN.lowerStruct.size();
        {
            const int updateColAlignment = origSN.size % childGridHeight;
            const int updateRowAlignment = origSN.size % childGridWidth;
            const int updateColShift = 
                elemental::Shift
                ( childGridRow, updateColAlignment, childGridHeight );
            const int updateRowShift = 
                elemental::Shift
                ( childGridCol, updateRowAlignment, childGridWidth );
            const int updateLocalHeight = 
                elemental::LocalLength
                ( updateSize, updateColShift, childGridHeight );
            const int updateLocalWidth = 
                elemental::LocalLength
                ( updateSize, updateRowShift, childGridWidth );
            for( int jChildLocal=0; 
                     jChildLocal<updateLocalWidth; ++jChildLocal )
            {
                const int jChild = updateRowShift + jChildLocal*childGridWidth;
                const int destGridCol = myChildRelIndices[jChild] % gridWidth;

                const int align = (jChild+updateRowAlignment) % childGridHeight;
                const int shift = 
                    (childGridRow+childGridHeight-align) % childGridHeight;
                const int localColShift = 
                    (jChild+shift-updateColShift) / childGridHeight;
                for( int iChildLocal=localColShift; 
                         iChildLocal<updateLocalHeight; ++iChildLocal )
                {
                    const int iChild = 
                        updateColShift + iChildLocal*childGridHeight;
                    const int destGridRow = 
                        myChildRelIndices[iChild] % gridHeight;

                    const int destRank = destGridRow + destGridCol*gridHeight;
                    ++factSN.numChildFactSendIndices[destRank];
                }
            }
        }

        // Fill numChildSolveSendIndices to use for many solves
        factSN.numChildSolveSendIndices.resize( teamSize );
        std::memset
        ( &factSN.numChildSolveSendIndices[0], 0, teamSize*sizeof(int) );
        {
            const int updateAlignment = origSN.size % childTeamSize;
            const int updateShift = 
                elemental::Shift
                ( childTeamRank, updateAlignment, childTeamSize );
            const int updateLocalHeight = 
                elemental::LocalLength
                ( updateSize, updateShift, childTeamSize );
            for( int iChildLocal=updateShift; 
                     iChildLocal<updateLocalHeight; ++iChildLocal )
            {
                const int iChild = updateShift + iChildLocal*childTeamSize;
                const int destRank = myChildRelIndices[iChild] % teamSize;
                ++factSN.numChildSolveSendIndices[destRank];
            }
        }

        // Fill {left,right}ChildFact{Col,Row}Indices so that we can reuse them
        // to compute our recv information for use in many factorizations
        factSN.leftChildFactColIndices.clear();
        for( int i=0; i<numLeftIndices; ++i )
            if( factSN.leftChildRelIndices[i] % gridHeight == gridRow )
                factSN.leftChildFactColIndices.push_back( i );
        factSN.leftChildFactRowIndices.clear();
        for( int i=0; i<numLeftIndices; ++i )
            if( factSN.leftChildRelIndices[i] % gridWidth == gridCol )
                factSN.leftChildFactRowIndices.push_back( i );
        factSN.rightChildFactColIndices.clear();
        for( int i=0; i<numRightIndices; ++i )
            if( factSN.rightChildRelIndices[i] % gridHeight == gridRow )
                factSN.rightChildFactColIndices.push_back( i );
        factSN.rightChildFactRowIndices.clear();
        for( int i=0; i<numRightIndices; ++i )
            if( factSN.rightChildRelIndices[i] % gridWidth == gridCol )
                factSN.rightChildFactRowIndices.push_back( i );

        // Fill {left,right}ChildSolveIndices for use in many solves
        factSN.leftChildSolveIndices.clear();
        for( int i=0; i<numLeftIndices; ++i )
            if( factSN.leftChildRelIndices[i] % teamSize == teamRank )
                factSN.leftChildSolveIndices.push_back( i );
        factSN.rightChildSolveIndices.clear();
        for( int i=0; i<numRightIndices; ++i )
            if( factSN.rightChildRelIndices[i] % teamSize == teamRank )
                factSN.rightChildSolveIndices.push_back( i );

        //
        // Compute the solve recv indices
        //
        const int leftChildTeamSize = teamSize / 2;
        const int rightChildTeamSize = teamSize / 2;
        factSN.childSolveRecvIndices.clear();
        factSN.childSolveRecvIndices.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftUpdateAlignment = 
            factSN.leftChildSize % leftChildTeamSize;
        for( int iPre=0; iPre<factSN.leftChildSolveIndices.size(); ++iPre )
        {
            const int iChild = factSN.leftChildSolveIndices[iPre];
            const int iFront = factSN.leftChildRelIndices[iChild];
            const int iFrontLocal = (iFront-teamRank) / teamSize;

            const int childRank = 
                (iChild+leftUpdateAlignment) % leftChildTeamSize;
            const int frontRank = childRank;
            factSN.childSolveRecvIndices[frontRank].push_back(iFrontLocal);
        }

        // Compute the recv indices for the right child
        const int rightUpdateAlignment = 
            factSN.rightChildSize % rightChildTeamSize;
        for( int iPre=0; iPre<factSN.rightChildSolveIndices.size(); ++iPre )
        {
            const int iChild = factSN.rightChildSolveIndices[iPre];
            const int iFront = factSN.rightChildRelIndices[iChild];
            const int iFrontLocal = (iFront-teamRank) / teamSize;

            const int childRank = 
                (iChild+rightUpdateAlignment) % rightChildTeamSize;
            const int frontRank = leftChildTeamSize + childRank;
            factSN.childSolveRecvIndices[frontRank].push_back(iFrontLocal);
        }

        // Optionally compute the recv indices for the factorization. 
        // This is optional since it requires a nontrivial amount of storage.
        if( storeFactRecvIndices )
            ComputeFactRecvIndices( factSN );
        else
            factSN.childFactRecvIndices.clear();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void clique::symbolic::ComputeFactRecvIndices( const DistSymmFactSupernode& sn )
{
#ifndef RELEASE
    PushCallStack("symbolic::ComputeFactRecvIndices");
#endif
    const int commRank = mpi::CommRank( sn.comm );
    const int commSize = mpi::CommSize( sn.comm );
    const int gridHeight = sn.gridHeight;
    const int gridWidth = commSize / gridHeight;
    const int gridRow = commRank % gridHeight;
    const int gridCol = commRank / gridHeight;

    // Compute the child grid dimensions (this could be stored in the supernode
    // if we generalize from powers of two).
    const int rightChildOffset = commSize / 2;
    int leftChildGridHeight, leftChildGridWidth,
        rightChildGridHeight, rightChildGridWidth;
    if( gridHeight >= gridWidth )
    {
        leftChildGridHeight = rightChildGridHeight = gridHeight / 2;
        leftChildGridWidth = rightChildGridWidth = gridWidth;
    }
    else
    {
        leftChildGridHeight = rightChildGridHeight = gridHeight;
        leftChildGridWidth = rightChildGridWidth = gridWidth / 2;
    }

    sn.childFactRecvIndices.clear();
    sn.childFactRecvIndices.resize( commSize );
    std::deque<int>::const_iterator it;

    // Compute the recv indices of the left child from each process 
    const int leftUpdateColAlignment = sn.leftChildSize % leftChildGridHeight;
    const int leftUpdateRowAlignment = sn.leftChildSize % leftChildGridWidth;
    for( int jPre=0; jPre<sn.leftChildFactRowIndices.size(); ++jPre )
    {
        const int jChild = sn.leftChildFactRowIndices[jPre];
        const int jFront = sn.leftChildRelIndices[jChild];
        const int jFrontLocal = (jFront-gridCol) / gridWidth;

        const int childCol = 
            (jChild+leftUpdateRowAlignment) % leftChildGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( sn.leftChildFactColIndices.begin(),
               sn.leftChildFactColIndices.end(), jChild );
        const int iPreStart = int(it-sn.leftChildFactColIndices.begin());
        for( int iPre=iPreStart; 
                 iPre<sn.leftChildFactColIndices.size(); ++iPre )
        {
            const int iChild = sn.leftChildFactColIndices[iPre];
            const int iFront = sn.leftChildRelIndices[iChild];
            const int iFrontLocal = (iFront-gridRow) / gridHeight;

            const int childRow = 
                (iChild+leftUpdateColAlignment) % leftChildGridHeight;
            const int childRank = childRow + childCol*leftChildGridHeight;

            const int frontRank = childRank;
            sn.childFactRecvIndices[frontRank].push_back(iFrontLocal);
            sn.childFactRecvIndices[frontRank].push_back(jFrontLocal);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    const int rightUpdateColAlignment = 
        sn.rightChildSize % rightChildGridHeight;
    const int rightUpdateRowAlignment = sn.rightChildSize % rightChildGridWidth;
    for( int jPre=0; jPre<sn.rightChildFactRowIndices.size(); ++jPre )
    {
        const int jChild = sn.rightChildFactRowIndices[jPre];
        const int jFront = sn.rightChildRelIndices[jChild];
        const int jFrontLocal = (jFront-gridCol) / gridWidth;

        const int childCol = 
            (jChild+rightUpdateRowAlignment) % rightChildGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( sn.rightChildFactColIndices.begin(),
               sn.rightChildFactColIndices.end(), jChild );
        const int iPreStart = int(it-sn.rightChildFactColIndices.begin());
        for( int iPre=iPreStart; 
                 iPre<sn.rightChildFactColIndices.size(); ++iPre )
        {
            const int iChild = sn.rightChildFactColIndices[iPre];
            const int iFront = sn.rightChildRelIndices[iChild];
            const int iFrontLocal = (iFront-gridRow) / gridHeight;

            const int childRow = 
                (iChild+rightUpdateColAlignment) % rightChildGridHeight;
            const int childRank = childRow + childCol*rightChildGridHeight;

            const int frontRank = rightChildOffset + childRank;
            sn.childFactRecvIndices[frontRank].push_back(iFrontLocal);
            sn.childFactRecvIndices[frontRank].push_back(jFrontLocal);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

