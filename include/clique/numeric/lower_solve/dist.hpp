/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERSOLVE_DIST_HPP
#define CLIQ_NUMERIC_LOWERSOLVE_DIST_HPP

namespace cliq {

template<typename F> 
void DistLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F> 
void DistLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

template<typename F>
void DistLowerBackwardSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F>
void DistLowerBackwardSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistLowerForwardSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    const bool frontsAre1d = FrontsAre1d( frontType );
    if( frontType != LDL_1D && 
        frontType != LDL_SELINV_1D && 
        frontType != LDL_SELINV_2D && 
        frontType != BLOCK_LDL_2D )
        throw std::logic_error("This solve mode is not yet implemented");

    // Copy the information from the local portion into the distributed leaf
    const SymmFront<F>& localRootFront = L.localFronts.back();
    const DistSymmFront<F>& distLeafFront = L.distFronts[0];
    const Grid& leafGrid = ( frontsAre1d ? distLeafFront.front1dL.Grid() 
                                         : distLeafFront.front2dL.Grid() );
    distLeafFront.work1d.LockedAttach
    ( localRootFront.work.Height(), localRootFront.work.Width(), 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim(), 
      leafGrid );
    
    // Perform the distributed portion of the forward solve
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = ( frontsAre1d ? front.front1dL.Grid()
                                         : front.front2dL.Grid() );
        mpi::Comm comm = grid.VCComm();
        const int commSize = mpi::CommSize( comm );

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmFront<F>& childFront = L.distFronts[s-1];
        const Grid& childGrid = ( frontsAre1d ? childFront.front1dL.Grid()
                                              : childFront.front2dL.Grid() );
        mpi::Comm childComm = childGrid.VCComm();
        const int childCommSize = mpi::CommSize( childComm );

        // Set up a workspace
        const int frontHeight = ( frontsAre1d ? front.front1dL.Height()
                                              : front.front2dL.Height() );
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( frontHeight, width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, node.size );
        WT = X.distNodes[s-1];
        elem::MakeZeros( WB );

        // Pack our child's update
        const MultiVecCommMeta& commMeta = node.multiVecMeta;
        DistMatrix<F,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childNode.size;
        DistMatrix<F,VC,STAR> childUpdate( childW.Grid() );
        LockedView( childUpdate, childW, childNode.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = commMeta.numChildSendInd[proc]*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const bool onLeft = childNode.onLeft;
        const std::vector<int>& myChildRelInd = 
            ( onLeft ? node.leftRelInd : node.rightRelInd );
        const int colShift = childUpdate.ColShift();
        const int localHeight = childUpdate.LocalHeight();
        std::vector<int> packOffsets = sendDispls;
        for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const int iChild = colShift + iChildLoc*childCommSize;
            const int destRank = myChildRelInd[iChild] % commSize;
            for( int jChild=0; jChild<width; ++jChild )
                sendBuffer[packOffsets[destRank]++] = 
                    childUpdate.GetLocal(iChildLoc,jChild);
        }
        std::vector<int>().swap( packOffsets );
        childW.Empty();
        if( s == 1 )
            L.localFronts.back().work.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int recvSize = commMeta.childRecvInd[proc].size()*width;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        VerifySendsAndRecvs( sendCounts, recvCounts, comm );
#endif

        // AllToAll to send and receive the child updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, comm );
        std::vector<F>().swap( sendBuffer );
        std::vector<int>().swap( sendCounts );
        std::vector<int>().swap( sendDispls );

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::vector<int>& recvInd = commMeta.childRecvInd[proc];
            for( unsigned k=0; k<recvInd.size(); ++k )
            {
                const int iFrontLoc = recvInd[k];
                const F* recvRow = &recvValues[k*width];
                F* WRow = W.Buffer( iFrontLoc, 0 );
                const int WLDim = W.LDim();
                for( int j=0; j<width; ++j )
                    WRow[j*WLDim] += recvRow[j];
            }
        }
        std::vector<F>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );

        // Now that the RHS is set up, perform this node's solve
        if( frontType == LDL_1D )
            FrontLowerForwardSolve( front.front1dL, W );
        else if( frontType == LDL_SELINV_1D )
            FrontFastLowerForwardSolve( front.front1dL, W );
        else if( frontType == LDL_SELINV_2D )
            FrontFastLowerForwardSolve( front.front2dL, W );
        else // frontType == BLOCK_LDL_2D
            FrontBlockLowerForwardSolve( front.front2dL, W );

        // Store this node's portion of the result
        X.distNodes[s-1] = WT;
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work1d.Empty();
}

template<typename F> 
inline void DistLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistLowerForwardSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    if( FrontsAre1d(frontType) )
        throw std::logic_error("1d solves not yet implemented");
    const bool computeCommMetas = ( X.commMetas.size() == 0 );
    if( computeCommMetas )
        X.ComputeCommMetas( info );

    // Copy the information from the local portion into the distributed leaf
    const SymmFront<F>& localRootFront = L.localFronts.back();
    const DistSymmFront<F>& distLeafFront = L.distFronts[0];
    const Grid& leafGrid = distLeafFront.front2dL.Grid();
    distLeafFront.work2d.LockedAttach
    ( localRootFront.work.Height(), localRootFront.work.Width(), 0, 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim(), 
      leafGrid );
    
    // Perform the distributed portion of the forward solve
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = front.front2dL.Grid();
        const int gridHeight = grid.Height();
        const int gridWidth = grid.Width();
        mpi::Comm comm = grid.VCComm();
        const int commSize = mpi::CommSize( comm );

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmFront<F>& childFront = L.distFronts[s-1];
        const Grid& childGrid = childFront.front2dL.Grid();
        const int childGridHeight = childGrid.Height();
        const int childGridWidth = childGrid.Width();

        // Set up a workspace
        const int frontHeight = front.front2dL.Height();
        DistMatrix<F>& W = front.work2d;
        W.SetGrid( grid );
        W.ResizeTo( frontHeight, width );
        DistMatrix<F> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, node.size );
        WT = X.distNodes[s-1];
        elem::MakeZeros( WB );

        // Pack our child's update
        const MatrixCommMeta& commMeta = X.commMetas[s-1];
        DistMatrix<F>& childW = childFront.work2d;
        const int updateSize = childW.Height()-childNode.size;
        DistMatrix<F> childUpdate( childW.Grid() );
        LockedView( childUpdate, childW, childNode.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = commMeta.numChildSendInd[proc];
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        // Pack send data
        const bool onLeft = childNode.onLeft;
        const std::vector<int>& myChildRelInd = 
            ( onLeft ? node.leftRelInd : node.rightRelInd );
        const int colShift = childUpdate.ColShift();
        const int rowShift = childUpdate.RowShift();
        const int localWidth = childUpdate.LocalWidth();
        const int localHeight = childUpdate.LocalHeight();
        std::vector<int> packOffsets = sendDispls;
        for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const int iChild = colShift + iChildLoc*childGridHeight;
            const int destRow = myChildRelInd[iChild] % gridHeight;
            for( int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
            {
                const int jChild = rowShift + jChildLoc*childGridWidth;
                const int destCol = jChild % gridWidth;
                const int destRank = destRow + destCol*gridHeight;
                sendBuffer[packOffsets[destRank]++] = 
                    childUpdate.GetLocal(iChildLoc,jChildLoc);
            }
        }
        std::vector<int>().swap( packOffsets );
        childW.Empty();
        if( s == 1 )
            L.localFronts.back().work.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int recvSize = commMeta.childRecvInd[proc].size()/2;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        VerifySendsAndRecvs( sendCounts, recvCounts, comm );
#endif

        // AllToAll to send and receive the child updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, comm );
        std::vector<F>().swap( sendBuffer );
        std::vector<int>().swap( sendCounts );
        std::vector<int>().swap( sendDispls );

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::vector<int>& recvInd = commMeta.childRecvInd[proc];
            for( unsigned k=0; k<recvInd.size()/2; ++k )
            {
                const int iFrontLoc = recvInd[2*k+0];
                const int jLoc = recvInd[2*k+1];
                W.UpdateLocal( iFrontLoc, jLoc, recvValues[k] );
            }
        }
        std::vector<F>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );

        // Now that the RHS is set up, perform this node's solve
        if( frontType == LDL_SELINV_2D )
            FrontFastLowerForwardSolve( front.front2dL, W );
        else if( frontType == LDL_2D )
            FrontLowerForwardSolve( front.front2dL, W );
        else // BLOCK_LDL_2D
            FrontBlockLowerForwardSolve( front.front2dL, W );

        // Store this node's portion of the result
        X.distNodes[s-1] = WT;
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work2d.Empty();
}

template<typename F>
inline void DistLowerBackwardSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistLowerBackwardSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    const bool frontsAre1d = FrontsAre1d( frontType );
    const bool blockLDL = ( frontType == BLOCK_LDL_2D );
    if( frontType != LDL_1D && 
        frontType != LDL_SELINV_1D && 
        frontType != LDL_SELINV_2D && 
        frontType != BLOCK_LDL_2D )
        throw std::logic_error("This solve mode is not yet implemented");

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    const SymmFront<F>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        View( localRootFront.work, X.localNodes.back() );
        if( !blockLDL )
            FrontLowerBackwardSolve
            ( orientation, localRootFront.frontL, localRootFront.work );
        else
            FrontBlockLowerBackwardSolve
            ( orientation, localRootFront.frontL, localRootFront.work );
    }
    else
    {
        const DistSymmFront<F>& rootFront = L.distFronts.back();
        View( rootFront.work1d, X.distNodes.back() );
        if( frontType == LDL_1D )
            FrontLowerBackwardSolve
            ( orientation, rootFront.front1dL, rootFront.work1d );
        else if( frontType == LDL_SELINV_1D )
            FrontFastLowerBackwardSolve
            ( orientation, rootFront.front1dL, rootFront.work1d );
        else if( frontType == LDL_SELINV_2D )
            FrontFastLowerBackwardSolve
            ( orientation, rootFront.front2dL, rootFront.work1d );
        else
            FrontBlockLowerBackwardSolve
            ( orientation, rootFront.front2dL, rootFront.work1d );
    }

    for( int s=numDistNodes-2; s>=0; --s )
    {
        const DistSymmNodeInfo& parentNode = info.distNodes[s+1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& parentFront = L.distFronts[s+1];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = ( frontsAre1d ? front.front1dL.Grid() 
                                         : front.front2dL.Grid() );
        const Grid& parentGrid = ( frontsAre1d ? parentFront.front1dL.Grid()
                                               : parentFront.front2dL.Grid() );
        mpi::Comm comm = grid.VCComm(); 
        mpi::Comm parentComm = parentGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int parentCommSize = mpi::CommSize( parentComm );
        const int frontHeight = ( frontsAre1d ? front.front1dL.Height()
                                              : front.front2dL.Height() );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( frontHeight, width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, node.size );
        Matrix<F>& XT = 
          ( s>0 ? X.distNodes[s-1].Matrix() : X.localNodes.back() );
        WT.Matrix() = XT;

        //
        // Set the bottom from the parent
        //

        // Pack the updates using the recv approach from the forward solve
        const MultiVecCommMeta& commMeta = parentNode.multiVecMeta;
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int sendSize = commMeta.childRecvInd[proc].size()*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendValues = &sendBuffer[sendDispls[proc]];
            const std::vector<int>& recvInd = commMeta.childRecvInd[proc];
            for( unsigned k=0; k<recvInd.size(); ++k )
            {
                const int iFrontLoc = recvInd[k];
                F* sendRow = &sendValues[k*width];
                const F* workRow = parentWork.LockedBuffer( iFrontLoc, 0 );
                const int workLDim = parentWork.LDim();
                for( int j=0; j<width; ++j )
                    sendRow[j] = workRow[j*workLDim];
            }
        }
        parentWork.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(parentCommSize), recvDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int recvSize = commMeta.numChildSendInd[proc]*width;
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        VerifySendsAndRecvs( sendCounts, recvCounts, parentComm );
#endif

        // AllToAll to send and recv parent updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, parentComm );
        std::vector<F>().swap( sendBuffer );
        std::vector<int>().swap( sendCounts );
        std::vector<int>().swap( sendDispls );

        // Unpack the updates using the send approach from the forward solve
        const bool onLeft = node.onLeft;
        const std::vector<int>& myRelInd = 
            ( onLeft ? parentNode.leftRelInd : parentNode.rightRelInd );
        const int colShift = WB.ColShift();
        const int localHeight = WB.LocalHeight();
        for( int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
        {
            const int iUpdate = colShift + iUpdateLoc*commSize;
            const int startRank = myRelInd[iUpdate] % parentCommSize;
            const F* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int j=0; j<width; ++j )
                WB.SetLocal(iUpdateLoc,j,recvBuf[j]);
            recvDispls[startRank] += width;
        }
        std::vector<F>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );

        // Call the custom node backward solve
        if( s > 0 )
        {
            if( frontType == LDL_1D )
                FrontLowerBackwardSolve
                ( orientation, front.front1dL, W );
            else if( frontType == LDL_SELINV_1D )
                FrontFastLowerBackwardSolve
                ( orientation, front.front1dL, W );
            else if( frontType == LDL_SELINV_2D )
                FrontFastLowerBackwardSolve
                ( orientation, front.front2dL, W );
            else // frontType == BLOCK_LDL_2D
                FrontBlockLowerBackwardSolve
                ( orientation, front.front2dL, front.work1d );
        }
        else
        {
            View( localRootFront.work, W.Matrix() );
            if( !blockLDL )
                FrontLowerBackwardSolve
                ( orientation, localRootFront.frontL, localRootFront.work );
            else
                FrontBlockLowerBackwardSolve
                ( orientation, localRootFront.frontL, localRootFront.work );
        }

        // Store this node's portion of the result
        XT = WT.Matrix();
    }
}

template<typename F>
inline void DistLowerBackwardSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistLowerBackwardSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    const SymmFrontType frontType = L.frontType;
    const bool blockLDL = ( frontType == BLOCK_LDL_2D );
    if( FrontsAre1d(frontType) )
        throw std::logic_error("1d solve mode is not yet implemented");

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    const SymmFront<F>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        View( localRootFront.work, X.localNodes.back() );
        if( !blockLDL )
            FrontLowerBackwardSolve
            ( orientation, localRootFront.frontL, localRootFront.work );
        else
            FrontBlockLowerBackwardSolve
            ( orientation, localRootFront.frontL, localRootFront.work );
    }
    else
    {
        const DistSymmFront<F>& rootFront = L.distFronts.back();
        View( rootFront.work2d, X.distNodes.back() );
        if( frontType == LDL_SELINV_2D )
            FrontFastLowerBackwardSolve
            ( orientation, rootFront.front2dL, rootFront.work2d );
        else if( frontType == LDL_2D )
            FrontLowerBackwardSolve
            ( orientation, rootFront.front2dL, rootFront.work2d );
        else
            FrontBlockLowerBackwardSolve
            ( orientation, rootFront.front2dL, rootFront.work2d );
    }

    for( int s=numDistNodes-2; s>=0; --s )
    {
        const DistSymmNodeInfo& parentNode = info.distNodes[s+1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& parentFront = L.distFronts[s+1];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& grid = front.front2dL.Grid();
        const int gridHeight = grid.Height(); 
        const int gridWidth = grid.Width();
        const Grid& parentGrid = parentFront.front2dL.Grid();
        const int parentGridHeight = parentGrid.Height();
        const int parentGridWidth = parentGrid.Width();
        mpi::Comm parentComm = parentGrid.VCComm();
        const int parentCommSize = mpi::CommSize( parentComm );
        const int frontHeight = front.front2dL.Height();

        // Set up a workspace
        DistMatrix<F>& W = front.work2d;
        W.SetGrid( grid );
        W.ResizeTo( frontHeight, width );
        DistMatrix<F> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, node.size );
        Matrix<F>& XT = 
          ( s>0 ? X.distNodes[s-1].Matrix() : X.localNodes.back() );
        WT.Matrix() = XT;

        //
        // Set the bottom from the parent
        //

        // Pack the updates using the recv approach from the forward solve
        const MatrixCommMeta& commMeta = X.commMetas[s];
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            // childRecvInd contains pairs of indices, but we will send one
            // floating-point value per pair
            const int sendSize = commMeta.childRecvInd[proc].size()/2;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F>& parentWork = parentFront.work2d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendValues = &sendBuffer[sendDispls[proc]];
            const std::vector<int>& recvInd = commMeta.childRecvInd[proc];
            for( unsigned k=0; k<recvInd.size()/2; ++k )
            {
                const int iFrontLoc = recvInd[2*k+0];
                const int jLoc = recvInd[2*k+1];
                sendValues[k] = parentWork.GetLocal( iFrontLoc, jLoc );
            }
        }
        parentWork.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(parentCommSize), recvDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int recvSize = commMeta.numChildSendInd[proc];
            recvCounts[proc] = recvSize;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += recvSize;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        VerifySendsAndRecvs( sendCounts, recvCounts, parentComm );
#endif

        // AllToAll to send and recv parent updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, parentComm );
        std::vector<F>().swap( sendBuffer );
        std::vector<int>().swap( sendCounts );
        std::vector<int>().swap( sendDispls );

        // Unpack the updates using the send approach from the forward solve
        const bool onLeft = node.onLeft;
        const std::vector<int>& myRelInd = 
            ( onLeft ? parentNode.leftRelInd : parentNode.rightRelInd );
        const int colShift = WB.ColShift();
        const int rowShift = WB.RowShift();
        const int localHeight = WB.LocalHeight();
        const int localWidth = WB.LocalWidth();
        for( int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
        {
            const int iUpdate = colShift + iUpdateLoc*gridHeight;
            const int startRow = myRelInd[iUpdate] % parentGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int jUpdate = rowShift + jLoc*gridWidth;
                const int startCol = jUpdate % parentGridWidth;
                const int startRank = startRow + startCol*parentGridHeight;
                WB.SetLocal
                ( iUpdateLoc, jLoc, recvBuffer[recvDispls[startRank]++] );
            }
        }
        std::vector<F>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );

        // Call the custom node backward solve
        if( s > 0 )
        {
            // TODO: Handle non-inverted case
            if( frontType == LDL_SELINV_2D )
                FrontFastLowerBackwardSolve
                ( orientation, front.front2dL, W );
            else if( frontType == LDL_2D )
                FrontLowerBackwardSolve
                ( orientation, front.front2dL, front.work2d );
            else
                FrontBlockLowerBackwardSolve
                ( orientation, front.front2dL, front.work2d );
        }
        else
        {
            View( localRootFront.work, W.Matrix() );
            if( !blockLDL )
                FrontLowerBackwardSolve
                ( orientation, localRootFront.frontL, localRootFront.work );
            else
                FrontBlockLowerBackwardSolve
                ( orientation, localRootFront.frontL, localRootFront.work );
        }

        // Store this node's portion of the result
        XT = WT.Matrix();
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERSOLVE_DIST_HPP
