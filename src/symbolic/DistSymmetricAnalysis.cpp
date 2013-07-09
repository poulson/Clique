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

inline void GetChildGridDims
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

inline void ComputeSolveMetadata1d
( const DistSymmElimTree& eTree, DistSymmInfo& info )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeSolveMetadata1d");
#endif
    // Handle the interface node
    info.distNodes[0].solveMeta1d.Empty();
    info.distNodes[0].solveMeta1d.localSize = info.localNodes.back().size;

    // Handle the truly distributed nodes
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNode& node = eTree.distNodes[s];
        DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        const int teamSize = mpi::CommSize( node.comm );
        const int teamRank = mpi::CommRank( node.comm );

        const DistSymmNode& childNode = eTree.distNodes[s-1];
        const DistSymmNodeInfo& childNodeInfo = info.distNodes[s-1];
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
            ( childNode.onLeft ? nodeInfo.leftRelInd : nodeInfo.rightRelInd );

        // Fill solveMeta1d.numChildSendInd
        SolveMetadata1d& solveMeta1d = nodeInfo.solveMeta1d;
        solveMeta1d.Empty();
        solveMeta1d.numChildSendInd.resize( teamSize );
        elem::MemZero( &solveMeta1d.numChildSendInd[0], teamSize );
        const int updateSize = childNodeInfo.lowerStruct.size();
        {
            const int align = childNodeInfo.size % childTeamSize;
            const int shift = Shift( childTeamRank, align, childTeamSize );
            const int localHeight = Length( updateSize, shift, childTeamSize );
            for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
            {
                const int iChild = shift + iChildLoc*childTeamSize;
                const int destRank = myRelInd[iChild] % teamSize;
                ++solveMeta1d.numChildSendInd[destRank];
            }
        }

        // Fill solveMeta1d.{left,right}Ind for use in many solves
        const int numLeftInd = nodeInfo.leftRelInd.size();
        const int numRightInd = nodeInfo.rightRelInd.size();
        for( int i=0; i<numLeftInd; ++i )
            if( nodeInfo.leftRelInd[i] % teamSize == teamRank )
                solveMeta1d.leftInd.push_back( i );
        for( int i=0; i<numRightInd; ++i )
            if( nodeInfo.rightRelInd[i] % teamSize == teamRank )
                solveMeta1d.rightInd.push_back( i );

        //
        // Compute the solve recv indices
        //
        solveMeta1d.childRecvInd.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftAlign = nodeInfo.leftSize % leftTeamSize;
        const int numLeftSolveInd = solveMeta1d.leftInd.size();
        for( int iPre=0; iPre<numLeftSolveInd; ++iPre )
        {
            const int iChild = solveMeta1d.leftInd[iPre];
            const int iFront = nodeInfo.leftRelInd[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+leftAlign) % leftTeamSize;
            const int frontRank = leftTeamOffset + childRank;
            solveMeta1d.childRecvInd[frontRank].push_back(iFrontLoc);
        }

        // Compute the recv indices for the right child
        const int rightAlign = nodeInfo.rightSize % rightTeamSize;
        const int numRightSolveInd = solveMeta1d.rightInd.size();
        for( int iPre=0; iPre<numRightSolveInd; ++iPre )
        {
            const int iChild = solveMeta1d.rightInd[iPre];
            const int iFront = nodeInfo.rightRelInd[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+rightAlign) % rightTeamSize;
            const int frontRank = rightTeamOffset + childRank;
            solveMeta1d.childRecvInd[frontRank].push_back(iFrontLoc);
        }

        solveMeta1d.localSize = Length(node.size,teamRank,teamSize);
    }
}

void ComputeSolveMetadata2d
( const DistSymmElimTree& eTree, DistSymmInfo& info, int width )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeSolveMetadata2d");
#endif
    // Handle the interface node
    info.distNodes[0].solveMeta2d.Empty();
    info.distNodes[0].solveMeta2d.localHeight = info.localNodes.back().size;
    info.distNodes[0].solveMeta2d.localWidth = width;

    // Handle the truly distributed nodes
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNode& node = eTree.distNodes[s];
        DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        const int teamSize = mpi::CommSize( node.comm );
        const int teamRank = mpi::CommRank( node.comm );
        const Grid& grid = *nodeInfo.grid;

        const DistSymmNode& childNode = eTree.distNodes[s-1];
        const DistSymmNodeInfo& childNodeInfo = info.distNodes[s-1];
        const int childTeamSize = mpi::CommSize( childNode.comm );
        const int childTeamRank = mpi::CommRank( childNode.comm );
        const bool inFirstTeam = ( childTeamRank == teamRank );
        const bool leftIsFirst = ( childNode.onLeft==inFirstTeam );
        const int leftTeamSize =
            ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
        const int rightTeamSize = teamSize - leftTeamSize;
        const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
        const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );
        const Grid& childGrid = *childNodeInfo.grid; 

        // Fill solveMeta2d.numChildSendInd
        SolveMetadata2d& solveMeta2d = nodeInfo.solveMeta2d;
        solveMeta2d.Empty();
        solveMeta2d.numChildSendInd.resize( teamSize );
        elem::MemZero( &solveMeta2d.numChildSendInd[0], teamSize );
        const int updateSize = childNodeInfo.lowerStruct.size();
        const std::vector<int>& myRelInd = 
            ( childNode.onLeft ? nodeInfo.leftRelInd : nodeInfo.rightRelInd );
        {
            const int colAlign = childNodeInfo.size % childGrid.Height();
            const int colShift = 
                Shift( childGrid.Row(), colAlign, childGrid.Height() );
            const int localHeight = 
                Length( updateSize, colShift, childGrid.Height() );

            const int rowAlign = 0;
            const int rowShift = 
                Shift( childGrid.Col(), rowAlign, childGrid.Width() );
            const int localWidth = Length( width, rowShift, childGrid.Width() );
            for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
            {
                const int iChild = colShift + iChildLoc*childGrid.Height();
                const int destRow = myRelInd[iChild] % grid.Height();
                for( int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
                {
                    const int jChild = rowShift + jChildLoc*childGrid.Width();
                    const int destCol = jChild % grid.Width();
                    const int destRank = destRow + destCol*grid.Height();
                    ++solveMeta2d.numChildSendInd[destRank];
                }
            }
        }

        // Fill solveMeta2d.{left,right}Ind for use in many solves
        const int numLeftInd = nodeInfo.leftRelInd.size();
        const int numRightInd = nodeInfo.rightRelInd.size();
        for( int i=0; i<numLeftInd; ++i )
            if( nodeInfo.leftRelInd[i] % grid.Height() == grid.Row() )
                solveMeta2d.leftRowInd.push_back( i );
        for( int i=0; i<numRightInd; ++i )
            if( nodeInfo.rightRelInd[i] % grid.Height() == grid.Row() )
                solveMeta2d.rightRowInd.push_back( i );

        // Get the child grid dimensions
        int childGridDims[4];
        GetChildGridDims( nodeInfo, childNodeInfo, childGridDims );
        const int leftGridHeight = childGridDims[0];
        const int leftGridWidth = childGridDims[1];
        const int rightGridHeight = childGridDims[2];
        const int rightGridWidth = childGridDims[3];

        //
        // Compute the solve recv indices
        //
        solveMeta2d.childRecvInd.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftColAlign = nodeInfo.leftSize % leftGridHeight;
        const int numLeftRowInd = solveMeta2d.leftRowInd.size();
        const int rowShift = Shift( grid.Col(), 0, grid.Width() );
        const int localWidth = Length( width, rowShift, grid.Width() );
        for( int iPre=0; iPre<numLeftRowInd; ++iPre )
        {
            const int iChild = solveMeta2d.leftRowInd[iPre];
            const int iFront = nodeInfo.leftRelInd[iChild];
            const int iFrontLoc = (iFront-grid.Row()) / grid.Height();
            const int childRow = (iChild+leftColAlign) % leftGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*grid.Width(); 
                const int childCol = j % leftGridWidth;
                const int childRank = childRow + childCol*leftGridHeight;
                const int frontRank = leftTeamOffset + childRank;
                solveMeta2d.childRecvInd[frontRank].push_back(iFrontLoc);
                solveMeta2d.childRecvInd[frontRank].push_back(jLoc);
             }
        }

        // Compute the recv indices for the right child
        const int rightColAlign = nodeInfo.rightSize % rightGridHeight;
        const int numRightRowInd = solveMeta2d.rightRowInd.size();
        for( int iPre=0; iPre<numRightRowInd; ++iPre )
        {
            const int iChild = solveMeta2d.rightRowInd[iPre];
            const int iFront = nodeInfo.rightRelInd[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;
            const int childRow = (iChild+rightColAlign) % rightGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*grid.Width();
                const int childCol = j % rightGridWidth;
                const int childRank = childRow + childCol*rightGridHeight;
                const int frontRank = rightTeamOffset + childRank;
                solveMeta2d.childRecvInd[frontRank].push_back(iFrontLoc);
                solveMeta2d.childRecvInd[frontRank].push_back(jLoc);
            }
        }

        solveMeta2d.localHeight = Length(node.size,grid.Row(),grid.Height());
        solveMeta2d.localWidth = localWidth;
    }
}

inline void ComputeFactorMetadata
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool computeFactRecvInd )
{
#ifndef RELEASE
    CallStackEntry entry("ComputeFactorMetadata");
#endif
    info.distNodes[0].factorMeta.Empty();
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNode& node = eTree.distNodes[s];
        DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        const int teamSize = mpi::CommSize( node.comm );

        const DistSymmNode& childNode = eTree.distNodes[s-1];
        const DistSymmNodeInfo& childNodeInfo = info.distNodes[s-1];

        // From ComputeStructAndRelInd
        FactorMetadata& factorMeta = nodeInfo.factorMeta;
        factorMeta.Empty();
        const int gridHeight = nodeInfo.grid->Height();
        const int gridWidth = nodeInfo.grid->Width();
        const int gridRow = nodeInfo.grid->Row();
        const int gridCol = nodeInfo.grid->Col();
        const int numLeftInd = nodeInfo.leftRelInd.size();
        const int numRightInd = nodeInfo.rightRelInd.size();
        for( int i=0; i<numLeftInd; ++i )
            if( nodeInfo.leftRelInd[i] % gridHeight == gridRow )
                factorMeta.leftColInd.push_back( i );
        for( int i=0; i<numLeftInd; ++i )
            if( nodeInfo.leftRelInd[i] % gridWidth == gridCol )
                factorMeta.leftRowInd.push_back( i );
        for( int i=0; i<numRightInd; ++i )
            if( nodeInfo.rightRelInd[i] % gridHeight == gridRow )
                factorMeta.rightColInd.push_back( i );
        for( int i=0; i<numRightInd; ++i )
            if( nodeInfo.rightRelInd[i] % gridWidth == gridCol )
                factorMeta.rightRowInd.push_back( i );

        // Fill factorMeta.numChildSendInd 
        const int childGridHeight = childNodeInfo.grid->Height();
        const int childGridWidth = childNodeInfo.grid->Width();
        const int childGridRow = childNodeInfo.grid->Row();
        const int childGridCol = childNodeInfo.grid->Col();
        const int mySize = childNodeInfo.size;
        const int updateSize = childNodeInfo.lowerStruct.size();
        factorMeta.numChildSendInd.resize( teamSize );
        elem::MemZero( &factorMeta.numChildSendInd[0], teamSize );
        const std::vector<int>& myRelInd = 
            ( childNode.onLeft ? nodeInfo.leftRelInd : nodeInfo.rightRelInd );
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
                    ++factorMeta.numChildSendInd[destRank];
                }
            }
        }

        // Optionally compute the recv indices for the factorization. 
        // This is optional since it requires a nontrivial amount of storage.
        if( computeFactRecvInd )
            ComputeFactRecvInd( nodeInfo, childNodeInfo );
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

    ComputeFactorMetadata( eTree, info, computeFactRecvInd );
    
    // This is thankfully independent of the number of right-hand sides,   
    // unlike 2d solve metadata
    ComputeSolveMetadata1d( eTree, info );
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

    const FactorMetadata& factorMeta = node.factorMeta;
    factorMeta.childRecvInd.resize( teamSize );
    const int numLeftColInd = factorMeta.leftColInd.size();
    const int numLeftRowInd = factorMeta.leftRowInd.size();
    const int numRightColInd = factorMeta.rightColInd.size();
    const int numRightRowInd = factorMeta.rightRowInd.size();
    std::deque<int>::const_iterator it;

    // Compute the recv indices of the left child from each process 
    const int gridRow = node.grid->Row();
    const int gridCol = node.grid->Col();
    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int childTeamRank = mpi::CommRank( childNode.comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
    for( int jPre=0; jPre<numLeftRowInd; ++jPre )
    {
        const int jChild = factorMeta.leftRowInd[jPre];
        const int jFront = node.leftRelInd[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid left jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.leftSize) % leftGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( factorMeta.leftColInd.begin(),
               factorMeta.leftColInd.end(), jChild );
        const int iPreStart = int(it-factorMeta.leftColInd.begin());

        for( int iPre=iPreStart; iPre<numLeftColInd; ++iPre )
        {
            const int iChild = factorMeta.leftColInd[iPre];
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
            factorMeta.childRecvInd[frontRank].push_back(iFrontLoc);
            factorMeta.childRecvInd[frontRank].push_back(jFrontLoc);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );
    for( int jPre=0; jPre<numRightRowInd; ++jPre )
    {
        const int jChild = factorMeta.rightRowInd[jPre];
        const int jFront = node.rightRelInd[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid right jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.rightSize) % rightGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( factorMeta.rightColInd.begin(),
               factorMeta.rightColInd.end(), jChild );
        const int iPreStart = int(it-factorMeta.rightColInd.begin());
        for( int iPre=iPreStart; iPre<numRightColInd; ++iPre )
        {
            const int iChild = factorMeta.rightColInd[iPre];
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
            factorMeta.childRecvInd[frontRank].push_back(iFrontLoc);
            factorMeta.childRecvInd[frontRank].push_back(jFrontLoc);
        }
    }
}

} // namespace cliq
