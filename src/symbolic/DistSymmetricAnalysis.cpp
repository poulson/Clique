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

inline void PairwiseExchangeLowerStruct
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

inline void BroadcastLowerStruct
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

inline void GetLowerStruct
( int& theirSize, std::vector<int>& theirLowerStruct,
  const DistSymmNode& node, const DistSymmNode& childNode, 
  const DistSymmNodeInfo& childNodeInfo )
{
    const int teamSize = mpi::CommSize( node.comm );
    const int childTeamSize = mpi::CommSize( childNode.comm );
    const int leftTeamSize =
        ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
    const int rightTeamSize = teamSize - leftTeamSize;
    if( leftTeamSize == rightTeamSize )
        PairwiseExchangeLowerStruct
        ( theirSize, theirLowerStruct, node, childNodeInfo );
    else
        BroadcastLowerStruct
        ( theirSize, theirLowerStruct, node, childNodeInfo );
}

inline void ComputeStructAndRelInd
( int theirSize, const std::vector<int>& theirLowerStruct,
  const DistSymmNode& node,         const DistSymmNode& childNode, 
        DistSymmNodeInfo& nodeInfo, const DistSymmNodeInfo& childNodeInfo )
{
    const std::vector<int>& myLowerStruct = childNodeInfo.lowerStruct;
#ifndef RELEASE
    if( !IsStrictlySorted(myLowerStruct) )
    {
        if( IsSorted(myLowerStruct) )
            throw std::logic_error("Repeat in my lower struct");
        else
            throw std::logic_error("My lower struct not sorted");
    }
    if( !IsStrictlySorted(theirLowerStruct) )
    {
        if( IsSorted(theirLowerStruct) )
            throw std::logic_error("Repeat in their lower struct");
        else
            throw std::logic_error("Their lower struct not sorted");
    }
    if( !IsStrictlySorted(node.lowerStruct) )
    {
        if( IsSorted(node.lowerStruct) )
            throw std::logic_error("Repeat in original struct");
        else
            throw std::logic_error("Original struct not sorted");
    }
#endif

    // Combine the children's structure
    std::vector<int> childrenStruct;
    Union( childrenStruct, myLowerStruct, theirLowerStruct );

    // Now add in the original lower structure
    std::vector<int> partialStruct;
    Union( partialStruct, childrenStruct, node.lowerStruct );

    // Now the node indices
    std::vector<int> nodeInd( node.size );
    for( int i=0; i<node.size; ++i )
        nodeInd[i] = node.offset + i;
    std::vector<int> fullStruct;
    Union( fullStruct, nodeInd, partialStruct );

    // Construct the relative indices of the original lower structure
    RelativeIndices
    ( nodeInfo.origLowerRelInd, node.lowerStruct, fullStruct );

    // Construct the relative indices of the children
    if( childNode.onLeft )
    {
        nodeInfo.leftSize = childNodeInfo.size;
        nodeInfo.rightSize = theirSize;
        RelativeIndices( nodeInfo.leftRelInd, myLowerStruct, fullStruct );
        RelativeIndices( nodeInfo.rightRelInd, theirLowerStruct, fullStruct );
    }
    else
    {
        nodeInfo.leftSize = theirSize;
        nodeInfo.rightSize = childNodeInfo.size;
        RelativeIndices( nodeInfo.leftRelInd, theirLowerStruct, fullStruct );
        RelativeIndices( nodeInfo.rightRelInd, myLowerStruct, fullStruct );
    }

    // Form lower structure of this node by removing the node indices
    const int teamRank = mpi::CommRank( node.comm );
    const int lowerStructSize = fullStruct.size() - node.size;
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
}

inline void ComputeMultiVecCommMeta( DistSymmInfo& info )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeMultiVecCommMeta");
#endif
    // Handle the interface node
    info.distNodes[0].multiVecMeta.Empty();
    info.distNodes[0].multiVecMeta.localSize = info.localNodes.back().size;

    // Handle the truly distributed nodes
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
    {
        DistSymmNodeInfo& node = info.distNodes[s];
        const int teamSize = mpi::CommSize( node.comm );
        const int teamRank = mpi::CommRank( node.comm );

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const int childTeamSize = mpi::CommSize( childNode.comm );
        const int childTeamRank = mpi::CommRank( childNode.comm );
        const bool inFirstTeam = ( childTeamRank == teamRank );
        const bool leftIsFirst = ( childNode.onLeft==inFirstTeam );
        const int leftTeamSize =
            ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
        const int rightTeamSize = teamSize - leftTeamSize;
        const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
        const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );

        const std::vector<int>& myRelInd = 
            ( childNode.onLeft ? node.leftRelInd : node.rightRelInd );

        // Fill numChildSendInd
        MultiVecCommMeta& commMeta = node.multiVecMeta;
        commMeta.Empty();
        commMeta.numChildSendInd.resize( teamSize );
        elem::MemZero( &commMeta.numChildSendInd[0], teamSize );
        const int updateSize = childNode.lowerStruct.size();
        {
            const int align = childNode.size % childTeamSize;
            const int shift = Shift( childTeamRank, align, childTeamSize );
            const int localHeight = Length( updateSize, shift, childTeamSize );
            for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
            {
                const int iChild = shift + iChildLoc*childTeamSize;
                const int destRank = myRelInd[iChild] % teamSize;
                ++commMeta.numChildSendInd[destRank];
            }
        }

        const int numLeftInd = node.leftRelInd.size();
        const int numRightInd = node.rightRelInd.size();
        std::vector<int> leftInd, rightInd; 
        for( int i=0; i<numLeftInd; ++i )
            if( node.leftRelInd[i] % teamSize == teamRank )
                leftInd.push_back( i );
        for( int i=0; i<numRightInd; ++i )
            if( node.rightRelInd[i] % teamSize == teamRank )
                rightInd.push_back( i );

        //
        // Compute the solve recv indices
        //
        commMeta.childRecvInd.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftAlign = node.leftSize % leftTeamSize;
        const int numLeftSolveInd = leftInd.size();
        for( int iPre=0; iPre<numLeftSolveInd; ++iPre )
        {
            const int iChild = leftInd[iPre];
            const int iFront = node.leftRelInd[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+leftAlign) % leftTeamSize;
            const int frontRank = leftTeamOffset + childRank;
            commMeta.childRecvInd[frontRank].push_back(iFrontLoc);
        }

        // Compute the recv indices for the right child
        const int rightAlign = node.rightSize % rightTeamSize;
        const int numRightSolveInd = rightInd.size();
        for( int iPre=0; iPre<numRightSolveInd; ++iPre )
        {
            const int iChild = rightInd[iPre];
            const int iFront = node.rightRelInd[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+rightAlign) % rightTeamSize;
            const int frontRank = rightTeamOffset + childRank;
            commMeta.childRecvInd[frontRank].push_back(iFrontLoc);
        }

        commMeta.localSize = Length(node.size,teamRank,teamSize);
    }
}

inline void ComputeFactorCommMeta( DistSymmInfo& info, bool computeFactRecvInd )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeFactorCommMeta");
#endif
    info.distNodes[0].factorMeta.Empty();
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
    {
        DistSymmNodeInfo& node = info.distNodes[s];
        const int teamSize = mpi::CommSize( node.comm );
        const DistSymmNodeInfo& childNode = info.distNodes[s-1];

        // Fill factorMeta.numChildSendInd 
        FactorCommMeta& commMeta = node.factorMeta;
        commMeta.Empty();
        const int gridHeight = node.grid->Height();
        const int gridWidth = node.grid->Width();
        const int childGridHeight = childNode.grid->Height();
        const int childGridWidth = childNode.grid->Width();
        const int childGridRow = childNode.grid->Row();
        const int childGridCol = childNode.grid->Col();
        const int mySize = childNode.size;
        const int updateSize = childNode.lowerStruct.size();
        commMeta.numChildSendInd.resize( teamSize );
        elem::MemZero( &commMeta.numChildSendInd[0], teamSize );
        const std::vector<int>& myRelInd = 
            ( childNode.onLeft ? node.leftRelInd : node.rightRelInd );
        {
            const int colAlign = mySize % childGridHeight;
            const int rowAlign = mySize % childGridWidth;
            const int colShift = 
                Shift( childGridRow, colAlign, childGridHeight );
            const int rowShift = 
                Shift( childGridCol, rowAlign, childGridWidth );
            const int localHeight = 
                Length( updateSize, colShift, childGridHeight );
            const int localWidth = 
                Length( updateSize, rowShift, childGridWidth );
            for( int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
            {
                const int jChild = rowShift + jChildLoc*childGridWidth;
                const int destGridCol = myRelInd[jChild] % gridWidth;

                int localColShift;
                if( colShift > jChild )
                    localColShift = 0;
                else if( (jChild-colShift) % childGridHeight == 0 )
                    localColShift = (jChild-colShift)/childGridHeight;
                else
                    localColShift = (jChild-colShift)/childGridHeight + 1;
                for( int iChildLoc=localColShift; 
                         iChildLoc<localHeight; ++iChildLoc )
                {
                    const int iChild = colShift + iChildLoc*childGridHeight;
                    const int destGridRow = myRelInd[iChild] % gridHeight;

                    const int destRank = destGridRow + destGridCol*gridHeight;
                    ++commMeta.numChildSendInd[destRank];
                }
            }
        }

        // Optionally compute the recv indices for the factorization. 
        // This is optional since it requires a nontrivial amount of storage.
        if( computeFactRecvInd )
            ComputeFactRecvInd( node, childNode );
    }
}

//
// This is the part of the analysis that requires fine-grain parallelism.
// For now, we will assume that the distributed part of the elimination 
// tree is binary.
//
void DistSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool computeFactRecvInd )
{
#ifndef RELEASE
    CallStackEntry entry("DistSymmetricAnalysis");
#endif
    const unsigned numDist = eTree.distNodes.size();
    info.distNodes.resize( numDist );

    // The bottom node was analyzed locally, so just copy its results over
    const SymmNodeInfo& topLocal = info.localNodes.back();
    DistSymmNodeInfo& bottomDist = info.distNodes[0];
    bottomDist.onLeft = eTree.distNodes[0].onLeft;
    mpi::CommDup( eTree.distNodes[0].comm, bottomDist.comm );
    bottomDist.grid = new Grid( bottomDist.comm );
    bottomDist.size = topLocal.size;
    bottomDist.offset = topLocal.offset;
    bottomDist.myOffset = topLocal.myOffset;
    bottomDist.lowerStruct = topLocal.lowerStruct;
    bottomDist.origLowerStruct = topLocal.origLowerStruct;
    bottomDist.origLowerRelInd = topLocal.origLowerRelInd;
    bottomDist.leftRelInd = topLocal.leftRelInd;
    bottomDist.rightRelInd = topLocal.rightRelInd;
    bottomDist.leftSize = -1; // not needed, could compute though
    bottomDist.rightSize = -1; // not needed, could compute though

    // Perform the distributed part of the symbolic factorization
    int myOffset = bottomDist.myOffset + bottomDist.size;
    for( unsigned s=1; s<numDist; ++s )
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
        GetLowerStruct
        ( theirSize, theirLowerStruct, node, childNode, childNodeInfo );

        // Perform one level of symbolic factorization and then compute
        // a wide variety of relative indices
        ComputeStructAndRelInd
        ( theirSize, theirLowerStruct, node, childNode, 
          nodeInfo, childNodeInfo );

        myOffset += nodeInfo.size;
    }

    ComputeFactorCommMeta( info, computeFactRecvInd );
    
    // This is thankfully independent of the number of right-hand sides,   
    // unlike the 2d equivalent
    ComputeMultiVecCommMeta( info );
}

// TODO: Simplify this implementation
void ComputeFactRecvInd
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeFactRecvInd");
#endif
    // Communicate to get the grid sizes
    int childGridDims[4];
    GetChildGridDims( node, childNode, childGridDims );
    const int leftGridHeight = childGridDims[0];
    const int leftGridWidth = childGridDims[1];
    const int rightGridHeight = childGridDims[2];
    const int rightGridWidth = childGridDims[3];

    const int teamSize = mpi::CommSize( node.comm );
    const int teamRank = mpi::CommRank( node.comm );
    const bool onLeft = childNode.onLeft;
    const int childTeamSize = mpi::CommSize( childNode.comm );
    const int leftTeamSize =
        ( onLeft ? childTeamSize : teamSize-childTeamSize );
    const int rightTeamSize = teamSize - leftTeamSize;
#ifndef RELEASE
    if( leftTeamSize != leftGridHeight*leftGridWidth )
        throw std::runtime_error("Computed left grid incorrectly");
    if( rightTeamSize != rightGridHeight*rightGridWidth )
        throw std::runtime_error("Computed right grid incorrectly");
#endif

    const FactorCommMeta& commMeta = node.factorMeta;
    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int gridRow = node.grid->Row();
    const int gridCol = node.grid->Col();
    const int numLeftInd = node.leftRelInd.size();
    const int numRightInd = node.rightRelInd.size();
    std::vector<int> leftRowInd, leftColInd, rightRowInd, rightColInd;
    for( int i=0; i<numLeftInd; ++i )
        if( node.leftRelInd[i] % gridHeight == gridRow )
            leftColInd.push_back( i );
    for( int i=0; i<numLeftInd; ++i )
        if( node.leftRelInd[i] % gridWidth == gridCol )
            leftRowInd.push_back( i );
    for( int i=0; i<numRightInd; ++i )
        if( node.rightRelInd[i] % gridHeight == gridRow )
            rightColInd.push_back( i );
    for( int i=0; i<numRightInd; ++i )
        if( node.rightRelInd[i] % gridWidth == gridCol )
            rightRowInd.push_back( i );

    // Compute the recv indices of the left child from each process 
    const int childTeamRank = mpi::CommRank( childNode.comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
    commMeta.childRecvInd.resize( teamSize );
    std::vector<int>::const_iterator it;
    const int numLeftColInd = leftColInd.size();
    const int numLeftRowInd = leftRowInd.size();
    const int numRightColInd = rightColInd.size();
    const int numRightRowInd = rightRowInd.size();
    for( int jPre=0; jPre<numLeftRowInd; ++jPre )
    {
        const int jChild = leftRowInd[jPre];
        const int jFront = node.leftRelInd[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid left jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.leftSize) % leftGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound( leftColInd.begin(), leftColInd.end(), jChild );
        const int iPreStart = int(it-leftColInd.begin());
        for( int iPre=iPreStart; iPre<numLeftColInd; ++iPre )
        {
            const int iChild = leftColInd[iPre];
            const int iFront = node.leftRelInd[iChild];
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
            commMeta.childRecvInd[frontRank].push_back(iFrontLoc);
            commMeta.childRecvInd[frontRank].push_back(jFrontLoc);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );
    for( int jPre=0; jPre<numRightRowInd; ++jPre )
    {
        const int jChild = rightRowInd[jPre];
        const int jFront = node.rightRelInd[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid right jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.rightSize) % rightGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound( rightColInd.begin(), rightColInd.end(), jChild );
        const int iPreStart = int(it-rightColInd.begin());
        for( int iPre=iPreStart; iPre<numRightColInd; ++iPre )
        {
            const int iChild = rightColInd[iPre];
            const int iFront = node.rightRelInd[iChild];
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
            commMeta.childRecvInd[frontRank].push_back(iFrontLoc);
            commMeta.childRecvInd[frontRank].push_back(jFrontLoc);
        }
    }
}

void GetChildGridDims
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode, 
  int* childGridDims )
{
    const bool onLeft = childNode.onLeft;
    const int childTeamRank = mpi::CommRank( childNode.comm );
    elem::MemZero( childGridDims, 4 );
    if( onLeft && childTeamRank == 0 )
    {
        childGridDims[0] = childNode.grid->Height();
        childGridDims[1] = childNode.grid->Width();
    }
    else if( !onLeft && childTeamRank == 0 )
    {
        childGridDims[2] = childNode.grid->Height();
        childGridDims[3] = childNode.grid->Width();
    }
    mpi::AllReduce( childGridDims, 4, mpi::SUM, node.comm );
}

} // namespace cliq
