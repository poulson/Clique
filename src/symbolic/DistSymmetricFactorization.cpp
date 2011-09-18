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
        bool storeRecvIndices )
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
    bottomDistSN.leftChildColIndices.clear();
    bottomDistSN.leftChildRowIndices.clear();
    bottomDistSN.rightChildColIndices.clear();
    bottomDistSN.rightChildRowIndices.clear();
    bottomDistSN.numChildSendIndices.clear();
    bottomDistSN.leftChildRecvIndices.clear();
    bottomDistSN.rightChildRecvIndices.clear();

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

        // SendRecv the message lengths
        const int myChildLowerStructSize = factChildSN.lowerStruct.size();
        int theirChildLowerStructSize;
        mpi::SendRecv
        ( &myChildLowerStructSize, 1, partner, 0,
          &theirChildLowerStructSize, 1, partner, 0, distOrig.comm );
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
            leftIndices = &sendBuffer[0];
            rightIndices = &recvBuffer[0];
            numLeftIndices = sendBuffer.size();
            numRightIndices = recvBuffer.size();
        }
        else
        {
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

        // Construct the send and recv information
        // HERE
        // TODO: Loop over our local matrix, incrementing the send counts,
        //       using the gridHeight and gridWidth variables

        // Form lower structure of this node by removing the supernode indices
        factSN.lowerStruct.resize( fullStructSize );
        it = std::set_difference
        ( fullStruct.begin(), fullStruct.end(),
          supernodeIndices.begin(), supernodeIndices.end(),
          factSN.lowerStruct.begin() );
        factSN.lowerStruct.resize( int(it-factSN.lowerStruct.begin()) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

