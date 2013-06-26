/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "clique.hpp"

namespace cliq {

namespace internal {

inline void GetLowerStructFromPartner
( int& theirSize, std::vector<int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childNodeInfo )
{
    // Determine our partner's rank for this exchange in node's communicator
    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const int childTeamRank = mpi::CommRank( childNodeInfo.comm );
    const int myTeamSize = mpi::CommSize( childNodeInfo.comm );
    const int otherTeamSize = teamSize - myTeamSize;
    const bool inFirstTeam = ( teamRank == childTeamRank );
    const int partner =
        ( inFirstTeam ? teamRank+myTeamSize : teamRank-otherTeamSize );

    // SendRecv the message lengths
    const int mySize = childNodeInfo.size;
    const int myLowerStructSize = childNodeInfo.lowerStruct.size();
    const int initialSends[2] = { mySize, myLowerStructSize };
    int initialRecvs[2];
    mpi::SendRecv
    ( initialSends, 2, partner, 0,
      initialRecvs, 2, partner, 0, node.comm );
    theirSize = initialRecvs[0];
    const int theirLowerStructSize = initialRecvs[1];

    // Perform the exchange
    theirLowerStruct.resize( theirLowerStructSize );
    mpi::SendRecv
    ( &childNodeInfo.lowerStruct[0], myLowerStructSize, partner, 0,
      &theirLowerStruct[0], theirLowerStructSize, partner, 0, node.comm );
}

inline void GetLowerStruct
( int& theirSize, std::vector<int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNodeInfo& childNodeInfo )
{
    // Determine our partner's rank for this exchange in node's communicator
    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const int childTeamRank = mpi::CommRank( childNodeInfo.comm );
    const int myTeamSize = mpi::CommSize( childNodeInfo.comm );
    const int otherTeamSize = teamSize - myTeamSize;
    const bool inFirstTeam = ( teamRank == childTeamRank );

    if( childTeamRank == 0 )
    {
        const int partner =
            ( inFirstTeam ? teamRank+myTeamSize : teamRank-otherTeamSize );

        // SendRecv the message lengths
        const int mySize = childNodeInfo.size;
        const int myLowerStructSize = childNodeInfo.lowerStruct.size();
        const int initialSends[2] = { mySize, myLowerStructSize };
        int initialRecvs[2];
        mpi::SendRecv
        ( initialSends, 2, partner, 0,
          initialRecvs, 2, partner, 0, node.comm );
        theirSize = initialRecvs[0];
        const int theirLowerStructSize = initialRecvs[1];

        // Perform the exchange
        theirLowerStruct.resize( theirLowerStructSize );
        mpi::SendRecv
        ( &childNodeInfo.lowerStruct[0], myLowerStructSize, partner, 0,
          &theirLowerStruct[0], theirLowerStructSize, partner, 0, node.comm );

        // Broadcast the other team's child's sizes
        mpi::Broadcast( initialRecvs, 2, 0, childNodeInfo.comm );

        // Broadcast the other team's child's lower struct
        mpi::Broadcast
        ( &theirLowerStruct[0], theirLowerStructSize, 0, childNodeInfo.comm );
    } 
    else
    {
        // Receive the other team's child's sizes
        int initialRecvs[2];
        mpi::Broadcast( initialRecvs, 2, 0, childNodeInfo.comm );
        theirSize = initialRecvs[0];
        const int theirLowerStructSize = initialRecvs[1];

        // Receive the other team's child's lower struct
        theirLowerStruct.resize( theirLowerStructSize );
        mpi::Broadcast
        ( &theirLowerStruct[0], theirLowerStructSize, 0, childNodeInfo.comm );
    }
}

inline void ComputeStructAndRelIndices
( int theirSize, const std::vector<int>& theirLowerStruct,
  const DistSymmNode& node,         const DistSymmNode& childNode, 
        DistSymmNodeInfo& nodeInfo, const DistSymmNodeInfo& childNodeInfo )
{
    const std::vector<int>& myLowerStruct = childNodeInfo.lowerStruct;
    const int mySize = childNodeInfo.size;
    const int myLowerStructSize = myLowerStruct.size();
    const int theirLowerStructSize = theirLowerStruct.size();

    std::vector<int>::iterator it;
    std::vector<int> childrenStruct, partialStruct, fullStruct, nodeIndices;

    // Union the two child lower structures
#ifndef RELEASE
    for( int i=1; i<myLowerStructSize; ++i )
    {
        const int thisIndex = myLowerStruct[i];
        const int lastIndex = myLowerStruct[i-1];
        if( thisIndex < lastIndex )
        {
            std::ostringstream msg;
            msg << "My child's lower struct was not sorted";
            throw std::logic_error( msg.str().c_str() );
        }
        else if( thisIndex == lastIndex )
        {
            std::ostringstream msg;
            msg << "My child's lower struct had repeated index, " << thisIndex;
            throw std::logic_error( msg.str().c_str() );
        }
    }
    for( int i=1; i<theirLowerStructSize; ++i )
    {
        const int thisIndex = theirLowerStruct[i];
        const int lastIndex = theirLowerStruct[i-1];
        if( thisIndex < lastIndex )
        {
            std::ostringstream msg;
            msg << "Their child's lower struct was not sorted";
            throw std::logic_error( msg.str().c_str() );
        }
        else if( thisIndex == lastIndex )
        {
            std::ostringstream msg;
            msg << "Their child's lower struct had repeated index, "
                << thisIndex;
            throw std::logic_error( msg.str().c_str() );
        }
    }
#endif
    childrenStruct.resize( myLowerStructSize+theirLowerStructSize );
    it = std::set_union
    ( myLowerStruct.begin(), myLowerStruct.end(),
      theirLowerStruct.begin(), theirLowerStruct.end(), 
      childrenStruct.begin() );
    const int childrenStructSize = int(it-childrenStruct.begin());
    childrenStruct.resize( childrenStructSize );

    // Union the lower structure of this node
#ifndef RELEASE
    for( unsigned i=1; i<node.lowerStruct.size(); ++i )
    {
        const int thisIndex = node.lowerStruct[i];
        const int lastIndex = node.lowerStruct[i-1];
        if( thisIndex < lastIndex )
        {
            std::ostringstream msg;
            msg << "Original struct was not sorted";
            throw std::logic_error( msg.str().c_str() );
        }
        else if( thisIndex == lastIndex )
        {
            std::ostringstream msg;
            msg << "Original struct had repeated index, " << thisIndex;
            throw std::logic_error( msg.str().c_str() );
        }
    }
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

    const int teamRank = mpi::CommRank( node.comm );
    const int teamSize = mpi::CommSize( node.comm );
    const bool onLeft = childNode.onLeft;

    // Construct the relative indices of the children
    int numLeftIndices, numRightIndices;
    const int *leftIndices, *rightIndices;
    if( onLeft )
    {
        nodeInfo.leftSize = mySize;
        nodeInfo.rightSize = theirSize;
        leftIndices = &myLowerStruct[0];
        rightIndices = &theirLowerStruct[0];
        numLeftIndices = myLowerStructSize;
        numRightIndices = theirLowerStructSize;
    }
    else
    {
        nodeInfo.leftSize = theirSize;
        nodeInfo.rightSize = mySize;
        leftIndices = &theirLowerStruct[0];
        rightIndices = &myLowerStruct[0];
        numLeftIndices = theirLowerStructSize;
        numRightIndices = myLowerStructSize;
    }
    nodeInfo.leftRelIndices.resize( numLeftIndices );
    it = fullStruct.begin();
    for( int i=0; i<numLeftIndices; ++i )
    {
        const int index = leftIndices[i];
        it = std::lower_bound( it, fullStruct.end(), index );
#ifndef RELEASE
        if( it == fullStruct.end() )
            throw std::logic_error("Relative index failed for left indices");
#endif
        nodeInfo.leftRelIndices[i] = int(it-fullStruct.begin());
    }
    nodeInfo.rightRelIndices.resize( numRightIndices );
    it = fullStruct.begin();
    for( int i=0; i<numRightIndices; ++i )
    {
        const int index = rightIndices[i];
        it = std::lower_bound( it, fullStruct.end(), index );
#ifndef RELEASE
        if( it == fullStruct.end() )
            throw std::logic_error("Relative index failed for right indices");
#endif
        nodeInfo.rightRelIndices[i] = int(it-fullStruct.begin());
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

    // Fill {left,right}Fact{Col,Row}Indices so that we can reuse them
    // to compute our recv information for use in many factorizations
    const unsigned gridHeight = nodeInfo.grid->Height();
    const unsigned gridWidth = nodeInfo.grid->Width();
    const unsigned gridRow = nodeInfo.grid->Row();
    const unsigned gridCol = nodeInfo.grid->Col();
    std::deque<int>().swap( nodeInfo.leftFactColIndices );
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftRelIndices[i] % gridHeight == gridRow )
            nodeInfo.leftFactColIndices.push_back( i );
    std::deque<int>().swap( nodeInfo.leftFactRowIndices );
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftRelIndices[i] % gridWidth == gridCol )
            nodeInfo.leftFactRowIndices.push_back( i );
    std::deque<int>().swap( nodeInfo.rightFactColIndices );
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightRelIndices[i] % gridHeight == gridRow )
            nodeInfo.rightFactColIndices.push_back( i );
    std::deque<int>().swap( nodeInfo.rightFactRowIndices );
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightRelIndices[i] % gridWidth == gridCol )
            nodeInfo.rightFactRowIndices.push_back( i );

    // Fill {left,right}SolveIndices for use in many solves
    std::deque<int>().swap( nodeInfo.leftSolveIndices );
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftRelIndices[i] % teamSize == teamRank )
            nodeInfo.leftSolveIndices.push_back( i );
    std::deque<int>().swap( nodeInfo.rightSolveIndices );
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightRelIndices[i] % teamSize == teamRank )
            nodeInfo.rightSolveIndices.push_back( i );
}

} // namespace internal

//
// This is the part of the analysis that requires fine-grain parallelism.
// For now, we will assume that the distributed part of the elimination 
// tree is binary.
//
void DistSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool storeFactRecvIndices )
{
#ifndef RELEASE
    CallStackEntry entry("DistSymmetricAnalysis");
#endif
    const unsigned numNodes = eTree.distNodes.size();
    info.distNodes.resize( numNodes );

    // The bottom node was analyzed locally, so just copy its results over
    const SymmNodeInfo& topLocal = info.localNodes.back();
    DistSymmNodeInfo& bottomDist = info.distNodes[0];
    bottomDist.onLeft = eTree.distNodes[0].onLeft;
    mpi::CommDup( eTree.distNodes[0].comm, bottomDist.comm );
    bottomDist.grid = new Grid( bottomDist.comm );
    bottomDist.size = topLocal.size;
    bottomDist.localSize1d = topLocal.size;
    bottomDist.offset = topLocal.offset;
    bottomDist.myOffset = topLocal.myOffset;
    bottomDist.localOffset1d = topLocal.myOffset;
    bottomDist.lowerStruct = topLocal.lowerStruct;
    bottomDist.origLowerStruct = topLocal.origLowerStruct;
    bottomDist.origLowerRelIndices = topLocal.origLowerRelIndices;
    bottomDist.leftRelIndices = topLocal.leftRelIndices;
    bottomDist.rightRelIndices = topLocal.rightRelIndices;
    bottomDist.leftSize = -1; // not needed, could compute though
    bottomDist.rightSize = -1; // not needed, could compute though
    std::deque<int>().swap( bottomDist.leftFactColIndices );
    std::deque<int>().swap( bottomDist.leftFactRowIndices );
    std::deque<int>().swap( bottomDist.rightFactColIndices );
    std::deque<int>().swap( bottomDist.rightFactRowIndices );
    std::vector<int>().swap( bottomDist.numChildFactSendIndices );
    std::vector<std::deque<int> >().swap( bottomDist.childFactRecvIndices );
    bottomDist.localOffset1d = topLocal.myOffset;
    std::deque<int>().swap( bottomDist.leftSolveIndices );
    std::deque<int>().swap( bottomDist.rightSolveIndices );
    std::vector<std::deque<int> >().swap( bottomDist.childSolveRecvIndices );

    // Perform the distributed part of the symbolic factorization
    int myOffset = bottomDist.myOffset + bottomDist.size;
    int localOffset1d = bottomDist.localOffset1d + bottomDist.size;
    for( unsigned s=1; s<numNodes; ++s )
    {
        const DistSymmNode& node = eTree.distNodes[s];
        const DistSymmNode& childNode = eTree.distNodes[s-1];
        const DistSymmNodeInfo& childNodeInfo = info.distNodes[s-1];
        DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        nodeInfo.onLeft = node.onLeft;
        nodeInfo.size = node.size;
        nodeInfo.offset = node.offset;
        nodeInfo.myOffset = myOffset;
        nodeInfo.origLowerStruct = node.lowerStruct;

        // Duplicate the communicator from the distributed eTree 
        mpi::CommDup( node.comm, nodeInfo.comm );
        nodeInfo.grid = new Grid( nodeInfo.comm );

        // Get the lower struct for the child we do not share
        int theirSize;
        std::vector<int> theirLowerStruct;
        const int teamSize = mpi::CommSize( node.comm );
        const int teamRank = mpi::CommRank( node.comm );
        const int childTeamRank = mpi::CommRank( childNode.comm );
        const int childTeamSize = mpi::CommSize( childNode.comm );
        const bool onLeft = childNode.onLeft;
        const int leftTeamSize = 
            ( onLeft ? childTeamSize : teamSize-childTeamSize );
        const int rightTeamSize = teamSize - leftTeamSize;
        const bool inFirstTeam = ( childTeamRank == teamRank );
        const bool leftIsFirst = ( onLeft==inFirstTeam );
        const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
        const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );
        if( leftTeamSize == rightTeamSize )
            internal::GetLowerStructFromPartner
            ( theirSize, theirLowerStruct, node, childNodeInfo );
        else
            internal::GetLowerStruct
            ( theirSize, theirLowerStruct, node, childNodeInfo );

        // Perform one level of symbolic factorization and then compute
        // a wide variety of relative indices
        internal::ComputeStructAndRelIndices
        ( theirSize, theirLowerStruct, node, childNode, 
          nodeInfo, childNodeInfo );

        // Fill numChildFactSendIndices so that we can reuse it for many facts.
        const unsigned gridHeight = nodeInfo.grid->Height();
        const unsigned gridWidth = nodeInfo.grid->Width();
        const unsigned childGridHeight = childNodeInfo.grid->Height();
        const unsigned childGridWidth = childNodeInfo.grid->Width();
        const unsigned childGridRow = childNodeInfo.grid->Row();
        const unsigned childGridCol = childNodeInfo.grid->Col();
        const int mySize = childNodeInfo.size;
        const int myLowerStructSize = childNodeInfo.lowerStruct.size();
        nodeInfo.numChildFactSendIndices.resize( teamSize );
        elem::MemZero( &nodeInfo.numChildFactSendIndices[0], teamSize );

        const std::vector<int>& myRelIndices = 
            ( onLeft ? nodeInfo.leftRelIndices : nodeInfo.rightRelIndices );
        const int updateSize = myLowerStructSize;
        {
            const int updateColAlignment = mySize % childGridHeight;
            const int updateRowAlignment = mySize % childGridWidth;
            const int updateColShift = 
                Shift<int>( childGridRow, updateColAlignment, childGridHeight );
            const int updateRowShift = 
                Shift<int>( childGridCol, updateRowAlignment, childGridWidth );
            const int updateLocalHeight = 
                Length<int>( updateSize, updateColShift, childGridHeight );
            const int updateLocalWidth = 
                Length<int>( updateSize, updateRowShift, childGridWidth );
            for( int jChildLoc=0; jChildLoc<updateLocalWidth; ++jChildLoc )
            {
                const int jChild = updateRowShift + jChildLoc*childGridWidth;
                const int destGridCol = myRelIndices[jChild] % gridWidth;

                int localColShift;
                if( updateColShift > jChild )
                    localColShift = 0;
                else if( (jChild-updateColShift) % childGridHeight == 0 )
                    localColShift = (jChild-updateColShift)/childGridHeight;
                else
                    localColShift = (jChild-updateColShift)/childGridHeight + 1;
                for( int iChildLoc=localColShift; 
                         iChildLoc<updateLocalHeight; ++iChildLoc )
                {
                    const int iChild = 
                        updateColShift + iChildLoc*childGridHeight;
                    const int destGridRow = myRelIndices[iChild] % gridHeight;

                    const int destRank = destGridRow + destGridCol*gridHeight;
                    ++nodeInfo.numChildFactSendIndices[destRank];
                }
            }
        }

        // Fill numChildSolveSendIndices to use for many solves
        nodeInfo.numChildSolveSendIndices.resize( teamSize );
        elem::MemZero( &nodeInfo.numChildSolveSendIndices[0], teamSize );
        {
            const int updateAlignment = mySize % childTeamSize;
            const int updateShift = 
                Shift<int>( childTeamRank, updateAlignment, childTeamSize );
            const int updateLocalHeight = 
                Length<int>( updateSize, updateShift, childTeamSize );
            for( int iChildLoc=0; 
                     iChildLoc<updateLocalHeight; ++iChildLoc )
            {
                const int iChild = updateShift + iChildLoc*childTeamSize;
                const int destRank = myRelIndices[iChild] % teamSize;
                ++nodeInfo.numChildSolveSendIndices[destRank];
            }
        }

        //
        // Compute the solve recv indices
        //
        std::vector<std::deque<int> >().swap( nodeInfo.childSolveRecvIndices );
        nodeInfo.childSolveRecvIndices.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftUpdateAlignment = nodeInfo.leftSize % leftTeamSize;
        const int numLeftSolveIndices = nodeInfo.leftSolveIndices.size();
        for( int iPre=0; iPre<numLeftSolveIndices; ++iPre )
        {
            const int iChild = nodeInfo.leftSolveIndices[iPre];
            const int iFront = nodeInfo.leftRelIndices[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+leftUpdateAlignment) % leftTeamSize;
            const int frontRank = leftTeamOffset + childRank;
            nodeInfo.childSolveRecvIndices[frontRank].push_back(iFrontLoc);
        }

        // Compute the recv indices for the right child
        const int rightUpdateAlignment = nodeInfo.rightSize % rightTeamSize;
        const int numRightSolveIndices = nodeInfo.rightSolveIndices.size();
        for( int iPre=0; iPre<numRightSolveIndices; ++iPre )
        {
            const int iChild = nodeInfo.rightSolveIndices[iPre];
            const int iFront = nodeInfo.rightRelIndices[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+rightUpdateAlignment) % rightTeamSize;
            const int frontRank = rightTeamOffset + childRank;
            nodeInfo.childSolveRecvIndices[frontRank].push_back(iFrontLoc);
        }

        // Optionally compute the recv indices for the factorization. 
        // This is optional since it requires a nontrivial amount of storage.
        if( storeFactRecvIndices )
            ComputeFactRecvIndices( nodeInfo, childNodeInfo );
        else
            std::vector<std::deque<int> >().swap
            ( nodeInfo.childFactRecvIndices );

        nodeInfo.localSize1d = Length<int>(node.size,teamRank,teamSize);
        nodeInfo.localOffset1d = localOffset1d;

        myOffset += nodeInfo.size;
        localOffset1d += nodeInfo.localSize1d;
    }
}

void ComputeFactRecvIndices
( const DistSymmNodeInfo& node, 
  const DistSymmNodeInfo& childNode )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeFactRecvIndices");
#endif
    const int teamSize = mpi::CommSize( node.comm );
    const int teamRank = mpi::CommRank( node.comm );
    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int gridRow = node.grid->Row();
    const int gridCol = node.grid->Col();
    const int childGridHeight = childNode.grid->Height();
    const int childGridWidth = childNode.grid->Width();
    const int childTeamRank = mpi::CommRank( childNode.comm );
    const int childTeamSize = mpi::CommSize( childNode.comm );
    const bool onLeft = childNode.onLeft;
    const int leftTeamSize =
        ( onLeft ? childTeamSize : teamSize-childTeamSize );
    const int rightTeamSize = teamSize - leftTeamSize;
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
    const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );

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
    if( leftTeamSize != leftGridHeight*leftGridWidth )
        throw std::runtime_error("Computed left grid incorrectly");
    if( rightTeamSize != rightGridHeight*rightGridWidth )
        throw std::runtime_error("Computed right grid incorrectly");
#endif

    std::vector<std::deque<int> >().swap( node.childFactRecvIndices );
    node.childFactRecvIndices.resize( teamSize );
    std::deque<int>::const_iterator it;

    // Compute the recv indices of the left child from each process 
    for( unsigned jPre=0; jPre<node.leftFactRowIndices.size(); ++jPre )
    {
        const int jChild = node.leftFactRowIndices[jPre];
        const int jFront = node.leftRelIndices[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid left jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.leftSize) % leftGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( node.leftFactColIndices.begin(),
               node.leftFactColIndices.end(), jChild );
        const int iPreStart = int(it-node.leftFactColIndices.begin());
        for( unsigned iPre=iPreStart; 
                      iPre<node.leftFactColIndices.size(); ++iPre )
        {
            const int iChild = node.leftFactColIndices[iPre];
            const int iFront = node.leftRelIndices[iChild];
#ifndef RELEASE
            if( iChild < jChild )
                throw std::logic_error("Invalid left iChild");
            if( (iFront-gridRow) % gridHeight != 0 )
                throw std::logic_error("Invalid left iFront");
#endif
            const int iFrontLoc = (iFront-gridRow) / gridHeight;

            const int childRow = (iChild+node.leftSize) % leftGridHeight;
            const int childRank = childRow + childCol*leftGridHeight;

            const int frontRank = leftTeamOffset + childRank;
            node.childFactRecvIndices[frontRank].push_back(iFrontLoc);
            node.childFactRecvIndices[frontRank].push_back(jFrontLoc);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    for( unsigned jPre=0; jPre<node.rightFactRowIndices.size(); ++jPre )
    {
        const int jChild = node.rightFactRowIndices[jPre];
        const int jFront = node.rightRelIndices[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid right jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.rightSize) % rightGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( node.rightFactColIndices.begin(),
               node.rightFactColIndices.end(), jChild );
        const int iPreStart = int(it-node.rightFactColIndices.begin());
        for( unsigned iPre=iPreStart; 
                      iPre<node.rightFactColIndices.size(); ++iPre )
        {
            const int iChild = node.rightFactColIndices[iPre];
            const int iFront = node.rightRelIndices[iChild];
#ifndef RELEASE
            if( iChild < jChild )
                throw std::logic_error("Invalid right iChild");
            if( (iFront-gridRow) % gridHeight != 0 )
                throw std::logic_error("Invalid right iFront");
#endif
            const int iFrontLoc = (iFront-gridRow) / gridHeight;

            const int childRow = (iChild+node.rightSize) % rightGridHeight;
            const int childRank = childRow + childCol*rightGridHeight;

            const int frontRank = rightTeamOffset + childRank;
            node.childFactRecvIndices[frontRank].push_back(iFrontLoc);
            node.childFactRecvIndices[frontRank].push_back(jFrontLoc);
        }
    }
}

} // namespace cliq
