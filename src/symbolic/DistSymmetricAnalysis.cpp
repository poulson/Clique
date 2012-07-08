/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
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

namespace cliq {

namespace internal {

inline void GetLowerStructFromPartner
( int leftChildTeamSize, 
  int& theirChildSize, std::vector<int>& theirChildLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childNodeInfo )
{
    // Determine our partner's rank for this exchange in node's 
    // communicator
    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const bool onLeft = ( teamRank < leftChildTeamSize );
    const int partner =
        ( onLeft ? teamRank+leftChildTeamSize
                 : teamRank-leftChildTeamSize );

    // SendRecv the message lengths
    const int myChildSize = childNodeInfo.size;
    const int myChildLowerStructSize = childNodeInfo.lowerStruct.size();
    const int initialSends[2] = { myChildSize, myChildLowerStructSize };
    int initialRecvs[2];
    mpi::SendRecv
    ( initialSends, 2, partner, 0,
      initialRecvs, 2, partner, 0, node.comm );
    theirChildSize = initialRecvs[0];
    const int theirChildLowerStructSize = initialRecvs[1];

    // Perform the exchange
    theirChildLowerStruct.resize( theirChildLowerStructSize );
    mpi::SendRecv
    ( &childNodeInfo.lowerStruct[0], myChildLowerStructSize, partner, 0,
      &theirChildLowerStruct[0], theirChildLowerStructSize, partner, 0, 
      node.comm );
}

inline void GetLowerStruct
( int leftChildTeamSize, 
  int& theirChildSize, std::vector<int>& theirChildLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childNodeInfo )
{
    // Determine our partner's rank for this exchange in node's 
    // communicator
    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const int childTeamRank = mpi::CommRank( childNodeInfo.comm );
    const bool onLeft = ( teamRank < leftChildTeamSize );

    if( childTeamRank == 0 )
    {
        const int partner = ( onLeft ? leftChildTeamSize : 0 );

        // SendRecv the message lengths
        const int myChildSize = childNodeInfo.size;
        const int myChildLowerStructSize = childNodeInfo.lowerStruct.size();
        const int initialSends[2] = { myChildSize, myChildLowerStructSize };
        int initialRecvs[2];
        mpi::SendRecv
        ( initialSends, 2, partner, 0,
          initialRecvs, 2, partner, 0, node.comm );
        theirChildSize = initialRecvs[0];
        const int theirChildLowerStructSize = initialRecvs[1];

        // Perform the exchange
        theirChildLowerStruct.resize( theirChildLowerStructSize );
        mpi::SendRecv
        ( &childNodeInfo.lowerStruct[0], myChildLowerStructSize, partner, 0,
          &theirChildLowerStruct[0], theirChildLowerStructSize, partner, 0, 
          node.comm );

        // Broadcast the other team's child's sizes
        mpi::Broadcast( initialRecvs, 2, 0, childNodeInfo.comm );

        // Broadcast the other team's child's lower struct
        mpi::Broadcast
        ( &theirChildLowerStruct[0], theirChildLowerStructSize, 
          0, childNodeInfo.comm );
    } 
    else
    {
        // Receive the other team's child's sizes
        int initialRecvs[2];
        mpi::Broadcast( initialRecvs, 2, 0, childNodeInfo.comm );
        theirChildSize = initialRecvs[0];
        const int theirChildLowerStructSize = initialRecvs[1];

        // Receive the other team's child's lower struct
        theirChildLowerStruct.resize( theirChildLowerStructSize );
        mpi::Broadcast
        ( &theirChildLowerStruct[0], theirChildLowerStructSize,
          0, childNodeInfo.comm );
    }
}

inline void ComputeStructAndRelIndices
( int leftChildTeamSize, 
  int theirChildSize, const std::vector<int>& theirChildLowerStruct,
  const DistSymmNode& node,         const DistSymmNode& childNode, 
        DistSymmNodeInfo& nodeInfo, const DistSymmNodeInfo& childNodeInfo )
{
    const std::vector<int>& myChildLowerStruct = childNodeInfo.lowerStruct;
    const int myChildSize = childNodeInfo.size;
    const int myChildLowerStructSize = myChildLowerStruct.size();
    const int theirChildLowerStructSize = theirChildLowerStruct.size();

    std::vector<int>::iterator it;
    std::vector<int> childrenStruct, partialStruct, fullStruct, nodeIndices;

    // Union the two child lower structures
#ifndef RELEASE
    for( int i=1; i<myChildLowerStructSize; ++i )
        if( myChildLowerStruct[i] <= myChildLowerStruct[i-1] )
            throw std::runtime_error("Child's lower struct was not sorted");
    for( int i=1; i<theirChildLowerStructSize; ++i )
        if( theirChildLowerStruct[i] <= theirChildLowerStruct[i-1] )
            throw std::runtime_error("Child's lower struct was not sorted");
#endif
    childrenStruct.resize
    ( myChildLowerStructSize+theirChildLowerStructSize );
    it = std::set_union
    ( myChildLowerStruct.begin(), myChildLowerStruct.end(),
      theirChildLowerStruct.begin(), theirChildLowerStruct.end(), 
      childrenStruct.begin() );
    const int childrenStructSize = int(it-childrenStruct.begin());
    childrenStruct.resize( childrenStructSize );

    // Union the lower structure of this node
#ifndef RELEASE
    for( unsigned i=1; i<node.lowerStruct.size(); ++i )
        if( node.lowerStruct[i] <= node.lowerStruct[i-1] )
            throw std::runtime_error("Original struct was not sorted");
#endif
    partialStruct.resize( childrenStructSize + node.lowerStruct.size() );
    it = std::set_union
    ( childrenStruct.begin(), childrenStruct.end(),
      node.lowerStruct.begin(), node.lowerStruct.end(),
      partialStruct.begin() );
    const int partialStructSize = int(it-partialStruct.begin());
    partialStruct.resize( partialStructSize );

    // Union again with the node indices
    nodeIndices.resize( node.size );
    for( int i=0; i<node.size; ++i )
        nodeIndices[i] = node.offset + i;
    fullStruct.resize( node.size + partialStructSize );
    it = std::set_union
    ( nodeIndices.begin(), nodeIndices.end(),
      partialStruct.begin(), partialStruct.end(),
      fullStruct.begin() );
    const int fullStructSize = int(it-fullStruct.begin());
    fullStruct.resize( fullStructSize );

    // Construct the relative indices of the original lower structure
    const int numOrigLowerIndices = node.lowerStruct.size();
    it = fullStruct.begin();
    nodeInfo.origLowerRelIndices.resize( numOrigLowerIndices );
    for( int i=0; i<numOrigLowerIndices; ++i )
    {
        const int index = node.lowerStruct[i];
        it = std::lower_bound( it, fullStruct.end(), index );
#ifndef RELEASE
        if( it == fullStruct.end() )
            throw std::logic_error("Relative index failed for original struct");
#endif
        nodeInfo.origLowerRelIndices[i] = int(it-fullStruct.begin());
    }

    // Construct the relative indices of the children
    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const bool onLeft = ( teamRank < leftChildTeamSize );
    int numLeftIndices, numRightIndices;
    const int *leftIndices, *rightIndices;
    if( onLeft )
    {
        nodeInfo.leftChildSize = myChildSize;
        nodeInfo.rightChildSize = theirChildSize;
        leftIndices = &myChildLowerStruct[0];
        rightIndices = &theirChildLowerStruct[0];
        numLeftIndices = myChildLowerStructSize;
        numRightIndices = theirChildLowerStructSize;
    }
    else
    {
        nodeInfo.leftChildSize = theirChildSize;
        nodeInfo.rightChildSize = myChildSize;
        leftIndices = &theirChildLowerStruct[0];
        rightIndices = &myChildLowerStruct[0];
        numLeftIndices = theirChildLowerStructSize;
        numRightIndices = myChildLowerStructSize;
    }
    nodeInfo.leftChildRelIndices.resize( numLeftIndices );
    it = fullStruct.begin();
    for( int i=0; i<numLeftIndices; ++i )
    {
        const int index = leftIndices[i];
        it = std::lower_bound( it, fullStruct.end(), index );
#ifndef RELEASE
        if( it == fullStruct.end() )
            throw std::logic_error("Relative index failed for left indices");
#endif
        nodeInfo.leftChildRelIndices[i] = int(it-fullStruct.begin());
    }
    nodeInfo.rightChildRelIndices.resize( numRightIndices );
    it = fullStruct.begin();
    for( int i=0; i<numRightIndices; ++i )
    {
        const int index = rightIndices[i];
        it = std::lower_bound( it, fullStruct.end(), index );
#ifndef RELEASE
        if( it == fullStruct.end() )
            throw std::logic_error("Relative index failed for right indices");
#endif
        nodeInfo.rightChildRelIndices[i] = int(it-fullStruct.begin());
    }

    // Form lower structure of this node by removing the node indices
    const int lowerStructSize = fullStructSize - node.size;
    nodeInfo.lowerStruct.resize( lowerStructSize );
    for( int i=0; i<lowerStructSize; ++i )
        nodeInfo.lowerStruct[i] = fullStruct[node.size+i];
#ifndef RELEASE
    // Ensure that the root process computed a lowerStruct of the same size
    int rootLowerStructSize;
    if( teamRank == 0 )
        rootLowerStructSize = lowerStructSize;
    mpi::Broadcast( &rootLowerStructSize, 1, 0, node.comm );
    if( rootLowerStructSize != lowerStructSize )
        throw std::runtime_error("Root has different lower struct size");
#endif

    // Fill {left,right}ChildFact{Col,Row}Indices so that we can reuse them
    // to compute our recv information for use in many factorizations
    const unsigned gridHeight = nodeInfo.grid->Height();
    const unsigned gridWidth = nodeInfo.grid->Width();
    const unsigned gridRow = nodeInfo.grid->Row();
    const unsigned gridCol = nodeInfo.grid->Col();
    nodeInfo.leftChildFactColIndices.clear();
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftChildRelIndices[i] % gridHeight == gridRow )
            nodeInfo.leftChildFactColIndices.push_back( i );
    nodeInfo.leftChildFactRowIndices.clear();
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftChildRelIndices[i] % gridWidth == gridCol )
            nodeInfo.leftChildFactRowIndices.push_back( i );
    nodeInfo.rightChildFactColIndices.clear();
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightChildRelIndices[i] % gridHeight == gridRow )
            nodeInfo.rightChildFactColIndices.push_back( i );
    nodeInfo.rightChildFactRowIndices.clear();
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightChildRelIndices[i] % gridWidth == gridCol )
            nodeInfo.rightChildFactRowIndices.push_back( i );

    // Fill {left,right}ChildSolveIndices for use in many solves
    nodeInfo.leftChildSolveIndices.clear();
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftChildRelIndices[i] % teamSize == teamRank )
            nodeInfo.leftChildSolveIndices.push_back( i );
    nodeInfo.rightChildSolveIndices.clear();
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightChildRelIndices[i] % teamSize == teamRank )
            nodeInfo.rightChildSolveIndices.push_back( i );
}

} // namespace internal

//
// This is the part of the analysis that requires fine-grain parallelism.
// For now, we will assume that the distributed part of the elimination 
// tree is binary.
//
void DistSymmetricAnalysis
( const SymmElimTree& eTree, SymmInfo& info, bool storeFactRecvIndices )
{
#ifndef RELEASE
    PushCallStack("DistSymmetricAnalysis");
#endif
    const unsigned numNodes = eTree.dist.nodes.size();
    info.dist.nodes.resize( numNodes );

    // The bottom node was analyzed locally, so just copy its results over
    const LocalSymmNodeInfo& topLocal = info.local.nodes.back();
    DistSymmNodeInfo& bottomDist = info.dist.nodes[0];
    mpi::CommDup( eTree.dist.nodes[0].comm, bottomDist.comm );
    bottomDist.grid = new Grid( bottomDist.comm );
    bottomDist.size = topLocal.size;
    bottomDist.localSize1d = topLocal.size;
    bottomDist.offset = topLocal.offset;
    bottomDist.myOffset = topLocal.myOffset;
    bottomDist.localOffset1d = topLocal.myOffset;
    bottomDist.lowerStruct = topLocal.lowerStruct;
    bottomDist.origLowerStruct = topLocal.origLowerStruct;
    bottomDist.origLowerRelIndices = topLocal.origLowerRelIndices;
    bottomDist.leftChildRelIndices = topLocal.leftChildRelIndices;
    bottomDist.rightChildRelIndices = topLocal.rightChildRelIndices;
    bottomDist.leftChildSize = -1; // not needed, could compute though
    bottomDist.rightChildSize = -1; // not needed, could compute though
    bottomDist.leftChildFactColIndices.clear();
    bottomDist.leftChildFactRowIndices.clear();
    bottomDist.rightChildFactColIndices.clear();
    bottomDist.rightChildFactRowIndices.clear();
    bottomDist.numChildFactSendIndices.clear();
    bottomDist.childFactRecvIndices.clear();
    bottomDist.localOffset1d = topLocal.myOffset;
    bottomDist.leftChildSolveIndices.clear();
    bottomDist.rightChildSolveIndices.clear();
    bottomDist.childSolveRecvIndices.clear();

    // Perform the distributed part of the symbolic factorization
    int myOffset = bottomDist.myOffset + bottomDist.size;
    int localOffset1d = bottomDist.localOffset1d + bottomDist.size;
    for( unsigned s=1; s<numNodes; ++s )
    {
        const DistSymmNode& node = eTree.dist.nodes[s];
        const DistSymmNode& childNode = eTree.dist.nodes[s-1];
        const DistSymmNodeInfo& childNodeInfo = info.dist.nodes[s-1];
        DistSymmNodeInfo& nodeInfo = info.dist.nodes[s];
        nodeInfo.size = node.size;
        nodeInfo.offset = node.offset;
        nodeInfo.myOffset = myOffset;
        nodeInfo.origLowerStruct = node.lowerStruct;

        // Duplicate the communicator from the distributed eTree 
        mpi::CommDup( node.comm, nodeInfo.comm );
        nodeInfo.grid = new Grid( nodeInfo.comm );

        // Get the lower struct for the child we do not share
        int theirChildSize;
        std::vector<int> theirChildLowerStruct;
        const int teamSize = mpi::CommSize( node.comm );
        const int teamRank = mpi::CommRank( node.comm );
        const int childTeamRank = mpi::CommRank( childNode.comm );
        const int childTeamSize = mpi::CommSize( childNode.comm );
        const bool onLeft = ( teamRank == childTeamRank );
        const int leftChildTeamSize = 
            ( onLeft ? childTeamSize : teamSize-childTeamSize );
        const int rightChildTeamSize = teamSize - leftChildTeamSize;
        if( leftChildTeamSize == rightChildTeamSize )
            internal::GetLowerStructFromPartner
            ( leftChildTeamSize, theirChildSize, theirChildLowerStruct, 
              node, childNodeInfo );
        else
            internal::GetLowerStruct
            ( leftChildTeamSize, theirChildSize, theirChildLowerStruct, 
              node, childNodeInfo );

        // Perform one level of symbolic factorization and then compute
        // a wide variety of relative indices
        internal::ComputeStructAndRelIndices
        ( leftChildTeamSize, theirChildSize, theirChildLowerStruct, 
          node, childNode, nodeInfo, childNodeInfo );

        // Fill numChildFactSendIndices so that we can reuse it for many facts.
        const unsigned gridHeight = nodeInfo.grid->Height();
        const unsigned gridWidth = nodeInfo.grid->Width();
        const unsigned childGridHeight = childNodeInfo.grid->Height();
        const unsigned childGridWidth = childNodeInfo.grid->Width();
        const unsigned childGridRow = childNodeInfo.grid->Row();
        const unsigned childGridCol = childNodeInfo.grid->Col();
        const int myChildSize = childNodeInfo.size;
        const int myChildLowerStructSize = childNodeInfo.lowerStruct.size();
        nodeInfo.numChildFactSendIndices.resize( teamSize );
        elem::MemZero( &nodeInfo.numChildFactSendIndices[0], teamSize );

        const std::vector<int>& myChildRelIndices = 
            ( onLeft ? nodeInfo.leftChildRelIndices 
                     : nodeInfo.rightChildRelIndices );
        const int updateSize = myChildLowerStructSize;
        {
            const int updateColAlignment = myChildSize % childGridHeight;
            const int updateRowAlignment = myChildSize % childGridWidth;
            const int updateColShift = 
                Shift<int>( childGridRow, updateColAlignment, childGridHeight );
            const int updateRowShift = 
                Shift<int>( childGridCol, updateRowAlignment, childGridWidth );
            const int updateLocalHeight = 
                LocalLength<int>( updateSize, updateColShift, childGridHeight );
            const int updateLocalWidth = 
                LocalLength<int>( updateSize, updateRowShift, childGridWidth );
            for( int jChildLocal=0; 
                     jChildLocal<updateLocalWidth; ++jChildLocal )
            {
                const int jChild = updateRowShift + jChildLocal*childGridWidth;
                const int destGridCol = myChildRelIndices[jChild] % gridWidth;

                int localColShift;
                if( updateColShift > jChild )
                    localColShift = 0;
                else if( (jChild-updateColShift) % childGridHeight == 0 )
                    localColShift = (jChild-updateColShift)/childGridHeight;
                else
                    localColShift = (jChild-updateColShift)/childGridHeight + 1;
                for( int iChildLocal=localColShift; 
                         iChildLocal<updateLocalHeight; ++iChildLocal )
                {
                    const int iChild = 
                        updateColShift + iChildLocal*childGridHeight;
                    const int destGridRow = 
                        myChildRelIndices[iChild] % gridHeight;

                    const int destRank = destGridRow + destGridCol*gridHeight;
                    ++nodeInfo.numChildFactSendIndices[destRank];
                }
            }
        }

        // Fill numChildSolveSendIndices to use for many solves
        nodeInfo.numChildSolveSendIndices.resize( teamSize );
        elem::MemZero( &nodeInfo.numChildSolveSendIndices[0], teamSize );
        {
            const int updateAlignment = myChildSize % childTeamSize;
            const int updateShift = 
                Shift<int>( childTeamRank, updateAlignment, childTeamSize );
            const int updateLocalHeight = 
                LocalLength<int>( updateSize, updateShift, childTeamSize );
            for( int iChildLocal=0; 
                     iChildLocal<updateLocalHeight; ++iChildLocal )
            {
                const int iChild = updateShift + iChildLocal*childTeamSize;
                const int destRank = myChildRelIndices[iChild] % teamSize;
                ++nodeInfo.numChildSolveSendIndices[destRank];
            }
        }

        //
        // Compute the solve recv indices
        //
        nodeInfo.childSolveRecvIndices.clear();
        nodeInfo.childSolveRecvIndices.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftUpdateAlignment = 
            nodeInfo.leftChildSize % leftChildTeamSize;
        const int numLeftSolveIndices = nodeInfo.leftChildSolveIndices.size();
        for( int iPre=0; iPre<numLeftSolveIndices; ++iPre )
        {
            const int iChild = nodeInfo.leftChildSolveIndices[iPre];
            const int iFront = nodeInfo.leftChildRelIndices[iChild];
            const int iFrontLocal = (iFront-teamRank) / teamSize;

            const int childRank = 
                (iChild+leftUpdateAlignment) % leftChildTeamSize;
            const int frontRank = childRank;
            nodeInfo.childSolveRecvIndices[frontRank].push_back(iFrontLocal);
        }

        // Compute the recv indices for the right child
        const int rightUpdateAlignment = 
            nodeInfo.rightChildSize % rightChildTeamSize;
        const int numRightSolveIndices = nodeInfo.rightChildSolveIndices.size();
        for( int iPre=0; iPre<numRightSolveIndices; ++iPre )
        {
            const int iChild = nodeInfo.rightChildSolveIndices[iPre];
            const int iFront = nodeInfo.rightChildRelIndices[iChild];
            const int iFrontLocal = (iFront-teamRank) / teamSize;

            const int childRank = 
                (iChild+rightUpdateAlignment) % rightChildTeamSize;
            const int frontRank = leftChildTeamSize + childRank;
            nodeInfo.childSolveRecvIndices[frontRank].push_back(iFrontLocal);
        }

        // Optionally compute the recv indices for the factorization. 
        // This is optional since it requires a nontrivial amount of storage.
        if( storeFactRecvIndices )
            ComputeFactRecvIndices( nodeInfo, childNodeInfo );
        else
            nodeInfo.childFactRecvIndices.clear();

        nodeInfo.localSize1d = LocalLength<int>(node.size,teamRank,teamSize);
        nodeInfo.localOffset1d = localOffset1d;

        myOffset += nodeInfo.size;
        localOffset1d += nodeInfo.localSize1d;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void ComputeFactRecvIndices
( const DistSymmNodeInfo& node, 
  const DistSymmNodeInfo& childNode )
{
#ifndef RELEASE
    PushCallStack("ComputeFactRecvIndices");
#endif
    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int gridRow = node.grid->Row();
    const int gridCol = node.grid->Col();
    const int childGridHeight = childNode.grid->Height();
    const int childGridWidth = childNode.grid->Width();
    const int childTeamRank = mpi::CommRank( childNode.comm );
    const int childTeamSize = mpi::CommSize( childNode.comm );
    const bool onLeft = ( teamRank == childTeamRank );
    const int leftChildTeamSize =
        ( onLeft ? childTeamSize : teamSize-childTeamSize );

    // Communicate to get the grid sizes
    int childGridDims[4] = { 0, 0, 0, 0 };
    if( onLeft && childTeamRank == 0 )
    {
        childGridDims[0] = childGridHeight;
        childGridDims[1] = childGridWidth;
    }
    else if( !onLeft && childTeamRank == 0 )
    {
        childGridDims[2] = childGridHeight;
        childGridDims[3] = childGridWidth;
    }
    mpi::AllReduce( childGridDims, 4, mpi::SUM, node.comm );
    const int leftGridHeight = childGridDims[0];
    const int leftGridWidth = childGridDims[1];
    const int rightGridHeight = childGridDims[2];
    const int rightGridWidth = childGridDims[3];

#ifndef RELEASE
    if( leftChildTeamSize != leftGridHeight*leftGridWidth )
    {
        std::cerr << "childGridDims: " 
                  << childGridDims[0] << " " << childGridDims[1] << " "
                  << childGridDims[2] << " " << childGridDims[3] << "\n"
                  << "leftChildTeamSize: " << leftChildTeamSize << std::endl;
        throw std::runtime_error("Computed left grid incorrectly");
    }
    const int rightChildTeamSize = teamSize - leftChildTeamSize;
    if( rightChildTeamSize != rightGridHeight*rightGridWidth )
        throw std::runtime_error("Computed right grid incorrectly");
#endif

    node.childFactRecvIndices.clear();
    node.childFactRecvIndices.resize( teamSize );
    std::deque<int>::const_iterator it;

    // Compute the recv indices of the left child from each process 
    for( unsigned jPre=0; jPre<node.leftChildFactRowIndices.size(); ++jPre )
    {
        const int jChild = node.leftChildFactRowIndices[jPre];
        const int jFront = node.leftChildRelIndices[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid left jFront");
#endif
        const int jFrontLocal = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.leftChildSize) % leftGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( node.leftChildFactColIndices.begin(),
               node.leftChildFactColIndices.end(), jChild );
        const int iPreStart = int(it-node.leftChildFactColIndices.begin());
        for( unsigned iPre=iPreStart; 
                      iPre<node.leftChildFactColIndices.size(); ++iPre )
        {
            const int iChild = node.leftChildFactColIndices[iPre];
            const int iFront = node.leftChildRelIndices[iChild];
#ifndef RELEASE
            if( iChild < jChild )
                throw std::logic_error("Invalid left iChild");
            if( (iFront-gridRow) % gridHeight != 0 )
                throw std::logic_error("Invalid left iFront");
#endif
            const int iFrontLocal = (iFront-gridRow) / gridHeight;

            const int childRow = (iChild+node.leftChildSize) % leftGridHeight;
            const int childRank = childRow + childCol*leftGridHeight;

            const int frontRank = childRank;
            node.childFactRecvIndices[frontRank].push_back(iFrontLocal);
            node.childFactRecvIndices[frontRank].push_back(jFrontLocal);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    for( unsigned jPre=0; jPre<node.rightChildFactRowIndices.size(); ++jPre )
    {
        const int jChild = node.rightChildFactRowIndices[jPre];
        const int jFront = node.rightChildRelIndices[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid right jFront");
#endif
        const int jFrontLocal = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.rightChildSize) % rightGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( node.rightChildFactColIndices.begin(),
               node.rightChildFactColIndices.end(), jChild );
        const int iPreStart = int(it-node.rightChildFactColIndices.begin());
        for( unsigned iPre=iPreStart; 
                 iPre<node.rightChildFactColIndices.size(); ++iPre )
        {
            const int iChild = node.rightChildFactColIndices[iPre];
            const int iFront = node.rightChildRelIndices[iChild];
#ifndef RELEASE
            if( iChild < jChild )
                throw std::logic_error("Invalid right iChild");
            if( (iFront-gridRow) % gridHeight != 0 )
                throw std::logic_error("Invalid right iFront");
#endif
            const int iFrontLocal = (iFront-gridRow) / gridHeight;

            const int childRow = (iChild+node.rightChildSize) % rightGridHeight;
            const int childRank = childRow + childCol*rightGridHeight;

            const int frontRank = leftChildTeamSize + childRank;
            node.childFactRecvIndices[frontRank].push_back(iFrontLocal);
            node.childFactRecvIndices[frontRank].push_back(jFrontLocal);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
