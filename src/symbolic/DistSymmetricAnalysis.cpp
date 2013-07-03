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

    // Union the two child lower structures
    std::vector<int> childrenStruct;
    Union( childrenStruct, myLowerStruct, theirLowerStruct );

    // Union the lower structure of this node
    std::vector<int> partialStruct;
    Union( partialStruct, childrenStruct, node.lowerStruct );

    // Union again with the node indices
    std::vector<int> nodeIndices( node.size );
    for( int i=0; i<node.size; ++i )
        nodeIndices[i] = node.offset + i;
    std::vector<int> fullStruct;
    Union( fullStruct, nodeIndices, partialStruct );

    // Construct the relative indices of the original lower structure
    RelativeIndices
    ( nodeInfo.origLowerRelIndices, node.lowerStruct, fullStruct );

    // Construct the relative indices of the children
    int numLeftIndices, numRightIndices;
    const int *leftIndices, *rightIndices;
    if( childNode.onLeft )
    {
        nodeInfo.leftSize = childNodeInfo.size;
        nodeInfo.rightSize = theirSize;
        leftIndices = &myLowerStruct[0];
        rightIndices = &theirLowerStruct[0];
        numLeftIndices = myLowerStruct.size();
        numRightIndices = theirLowerStruct.size();
        RelativeIndices
        ( nodeInfo.leftRelIndices, myLowerStruct, fullStruct );
        RelativeIndices
        ( nodeInfo.rightRelIndices, theirLowerStruct, fullStruct );
    }
    else
    {
        nodeInfo.leftSize = theirSize;
        nodeInfo.rightSize = childNodeInfo.size;
        leftIndices = &theirLowerStruct[0];
        rightIndices = &myLowerStruct[0];
        numLeftIndices = theirLowerStruct.size();
        numRightIndices = myLowerStruct.size();
        RelativeIndices
        ( nodeInfo.leftRelIndices, theirLowerStruct, fullStruct );
        RelativeIndices
        ( nodeInfo.rightRelIndices, myLowerStruct, fullStruct );
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

    // Fill factorMeta.{left,right}{Col,Row}Indices so that we can reuse them
    // to compute our recv information for use in many factorizations
    FactorMetadata& factorMeta = nodeInfo.factorMeta;
    const int gridHeight = nodeInfo.grid->Height();
    const int gridWidth = nodeInfo.grid->Width();
    const int gridRow = nodeInfo.grid->Row();
    const int gridCol = nodeInfo.grid->Col();
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftRelIndices[i] % gridHeight == gridRow )
            factorMeta.leftColIndices.push_back( i );
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftRelIndices[i] % gridWidth == gridCol )
            factorMeta.leftRowIndices.push_back( i );
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightRelIndices[i] % gridHeight == gridRow )
            factorMeta.rightColIndices.push_back( i );
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightRelIndices[i] % gridWidth == gridCol )
            factorMeta.rightRowIndices.push_back( i );

    // Fill solveMeta1d.{left,right}Indices for use in many solves
    SolveMetadata1d& solveMeta1d = nodeInfo.solveMeta1d;
    const int teamSize = mpi::CommSize( node.comm );
    for( int i=0; i<numLeftIndices; ++i )
        if( nodeInfo.leftRelIndices[i] % teamSize == teamRank )
            solveMeta1d.leftIndices.push_back( i );
    for( int i=0; i<numRightIndices; ++i )
        if( nodeInfo.rightRelIndices[i] % teamSize == teamRank )
            solveMeta1d.rightIndices.push_back( i );
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
    bottomDist.solveMeta1d.localSize = topLocal.size;
    bottomDist.offset = topLocal.offset;
    bottomDist.myOffset = topLocal.myOffset;
    bottomDist.lowerStruct = topLocal.lowerStruct;
    bottomDist.origLowerStruct = topLocal.origLowerStruct;
    bottomDist.origLowerRelIndices = topLocal.origLowerRelIndices;
    bottomDist.leftRelIndices = topLocal.leftRelIndices;
    bottomDist.rightRelIndices = topLocal.rightRelIndices;
    bottomDist.leftSize = -1; // not needed, could compute though
    bottomDist.rightSize = -1; // not needed, could compute though
    bottomDist.factorMeta.Empty();
    bottomDist.solveMeta1d.Empty();

    // Perform the distributed part of the symbolic factorization
    int myOffset = bottomDist.myOffset + bottomDist.size;
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
        //
        // TODO: Compute all metadata here?
        nodeInfo.factorMeta.Empty();
        nodeInfo.solveMeta1d.Empty();
        nodeInfo.solveMeta2d.Empty();
        internal::ComputeStructAndRelIndices
        ( theirSize, theirLowerStruct, node, childNode, 
          nodeInfo, childNodeInfo );

        // Fill factorMeta.numChildSendIndices 
        FactorMetadata& factorMeta = nodeInfo.factorMeta;
        const int gridHeight = nodeInfo.grid->Height();
        const int gridWidth = nodeInfo.grid->Width();
        const int childGridHeight = childNodeInfo.grid->Height();
        const int childGridWidth = childNodeInfo.grid->Width();
        const int childGridRow = childNodeInfo.grid->Row();
        const int childGridCol = childNodeInfo.grid->Col();
        const int mySize = childNodeInfo.size;
        const int myLowerStructSize = childNodeInfo.lowerStruct.size();
        factorMeta.numChildSendIndices.resize( teamSize );
        elem::MemZero( &factorMeta.numChildSendIndices[0], teamSize );
        const std::vector<int>& myRelIndices = 
            ( onLeft ? nodeInfo.leftRelIndices : nodeInfo.rightRelIndices );
        const int updateSize = myLowerStructSize;
        {
            const int updateColAlignment = mySize % childGridHeight;
            const int updateRowAlignment = mySize % childGridWidth;
            const int updateColShift = 
                Shift( childGridRow, updateColAlignment, childGridHeight );
            const int updateRowShift = 
                Shift( childGridCol, updateRowAlignment, childGridWidth );
            const int updateLocalHeight = 
                Length( updateSize, updateColShift, childGridHeight );
            const int updateLocalWidth = 
                Length( updateSize, updateRowShift, childGridWidth );
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
                    ++factorMeta.numChildSendIndices[destRank];
                }
            }
        }

        // Fill solveMeta1d.numChildSendIndices
        SolveMetadata1d& solveMeta1d = nodeInfo.solveMeta1d;
        solveMeta1d.numChildSendIndices.resize( teamSize );
        elem::MemZero( &solveMeta1d.numChildSendIndices[0], teamSize );
        {
            const int updateAlignment = mySize % childTeamSize;
            const int updateShift = 
                Shift( childTeamRank, updateAlignment, childTeamSize );
            const int updateLocalHeight = 
                Length( updateSize, updateShift, childTeamSize );
            for( int iChildLoc=0; 
                     iChildLoc<updateLocalHeight; ++iChildLoc )
            {
                const int iChild = updateShift + iChildLoc*childTeamSize;
                const int destRank = myRelIndices[iChild] % teamSize;
                ++solveMeta1d.numChildSendIndices[destRank];
            }
        }

        //
        // Compute the solve recv indices
        //
        solveMeta1d.childRecvIndices.resize( teamSize );

        // Compute the recv indices for the left child 
        const int leftUpdateAlignment = nodeInfo.leftSize % leftTeamSize;
        const int numLeftSolveIndices = solveMeta1d.leftIndices.size();
        for( int iPre=0; iPre<numLeftSolveIndices; ++iPre )
        {
            const int iChild = solveMeta1d.leftIndices[iPre];
            const int iFront = nodeInfo.leftRelIndices[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+leftUpdateAlignment) % leftTeamSize;
            const int frontRank = leftTeamOffset + childRank;
            solveMeta1d.childRecvIndices[frontRank].push_back(iFrontLoc);
        }

        // Compute the recv indices for the right child
        const int rightUpdateAlignment = nodeInfo.rightSize % rightTeamSize;
        const int numRightSolveIndices = solveMeta1d.rightIndices.size();
        for( int iPre=0; iPre<numRightSolveIndices; ++iPre )
        {
            const int iChild = solveMeta1d.rightIndices[iPre];
            const int iFront = nodeInfo.rightRelIndices[iChild];
            const int iFrontLoc = (iFront-teamRank) / teamSize;

            const int childRank = (iChild+rightUpdateAlignment) % rightTeamSize;
            const int frontRank = rightTeamOffset + childRank;
            solveMeta1d.childRecvIndices[frontRank].push_back(iFrontLoc);
        }

        // Optionally compute the recv indices for the factorization. 
        // This is optional since it requires a nontrivial amount of storage.
        if( storeFactRecvIndices )
            ComputeFactRecvIndices( nodeInfo, childNodeInfo );

        myOffset += nodeInfo.size;
        solveMeta1d.localSize = Length(node.size,teamRank,teamSize);
    }
}

void ComputeFactRecvIndices
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode )
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

    const FactorMetadata& factorMeta = node.factorMeta;
    factorMeta.childRecvIndices.resize( teamSize );
    std::deque<int>::const_iterator it;

    // Compute the recv indices of the left child from each process 
    for( unsigned jPre=0; jPre<factorMeta.leftRowIndices.size(); ++jPre )
    {
        const int jChild = factorMeta.leftRowIndices[jPre];
        const int jFront = node.leftRelIndices[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid left jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.leftSize) % leftGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( factorMeta.leftColIndices.begin(),
               factorMeta.leftColIndices.end(), jChild );
        const int iPreStart = int(it-factorMeta.leftColIndices.begin());
        for( unsigned iPre=iPreStart; 
                      iPre<factorMeta.leftColIndices.size(); ++iPre )
        {
            const int iChild = factorMeta.leftColIndices[iPre];
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
            factorMeta.childRecvIndices[frontRank].push_back(iFrontLoc);
            factorMeta.childRecvIndices[frontRank].push_back(jFrontLoc);
        }
    }
    
    // Compute the recv indices of the right child from each process 
    for( unsigned jPre=0; jPre<factorMeta.rightRowIndices.size(); ++jPre )
    {
        const int jChild = factorMeta.rightRowIndices[jPre];
        const int jFront = node.rightRelIndices[jChild];
#ifndef RELEASE
        if( (jFront-gridCol) % gridWidth != 0 )
            throw std::logic_error("Invalid right jFront");
#endif
        const int jFrontLoc = (jFront-gridCol) / gridWidth;
        const int childCol = (jChild+node.rightSize) % rightGridWidth;

        // Find the first iPre that maps to the lower triangle
        it = std::lower_bound
             ( factorMeta.rightColIndices.begin(),
               factorMeta.rightColIndices.end(), jChild );
        const int iPreStart = int(it-factorMeta.rightColIndices.begin());
        for( unsigned iPre=iPreStart; 
                      iPre<factorMeta.rightColIndices.size(); ++iPre )
        {
            const int iChild = factorMeta.rightColIndices[iPre];
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
            factorMeta.childRecvIndices[frontRank].push_back(iFrontLoc);
            factorMeta.childRecvIndices[frontRank].push_back(jFrontLoc);
        }
    }
}

} // namespace cliq
