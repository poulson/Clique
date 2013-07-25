/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERMULTIPLY_DIST_HPP
#define CLIQ_NUMERIC_LOWERMULTIPLY_DIST_HPP

namespace cliq {

template<typename T> 
void DistLowerMultiplyNormal
( int diagOffset, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X );

template<typename T> 
void DistLowerMultiplyTranspose
( Orientation orientation, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, 
  DistNodalMultiVec<T>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T> 
inline void DistLowerMultiplyNormal
( int diagOffset, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistLowerMultiplyNormal");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    if( L.frontType != SYMM_1D && L.frontType != LDL_1D )
        throw std::logic_error("This multiply mode is not yet implemented");

    // Copy the information from the local portion into the distributed leaf
    const SymmFront<T>& localRootFront = L.localFronts.back();
    const DistSymmFront<T>& distLeafFront = L.distFronts[0];
    LockedView
    ( distLeafFront.work1d,
      localRootFront.work.Height(), localRootFront.work.Width(), 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim(),
      distLeafFront.front1dL.Grid() );
    
    // Perform the distributed portion of the forward multiply
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<T>& childFront = L.distFronts[s-1];
        const DistSymmFront<T>& front = L.distFronts[s];
        const Grid& childGrid = childFront.front1dL.Grid();
        const Grid& grid = front.front1dL.Grid();
        mpi::Comm comm = grid.VCComm();
        mpi::Comm childComm = childGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int childCommSize = mpi::CommSize( childComm );

        // Set up a workspace
        DistMatrix<T,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( front.front1dL.Height(), width );
        DistMatrix<T,VC,STAR> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, node.size );
        WT = X.distNodes[s-1];
        elem::MakeZeros( WB );

        // Now that the right-hand side is set up, perform the multiply
        FrontLowerMultiply( NORMAL, diagOffset, front.front1dL, W );

        // Pack our child's update
        const MultiVecCommMeta& commMeta = node.multiVecMeta;
        DistMatrix<T,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childNode.size;
        DistMatrix<T,VC,STAR> childUpdate;
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
        std::vector<T> sendBuffer( sendBufferSize );

        const std::vector<int>& myChildRelInd = 
            ( childNode.onLeft ? node.leftRelInd : node.rightRelInd );
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
        std::vector<T> recvBuffer( recvBufferSize );
#ifndef RELEASE
        VerifySendsAndRecvs( sendCounts, recvCounts, comm );
#endif

        // AllToAll to send and receive the child updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, comm );
        std::vector<T>().swap( sendBuffer );
        std::vector<int>().swap( sendCounts );
        std::vector<int>().swap( sendDispls );

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const T* recvValues = &recvBuffer[recvDispls[proc]];
            const std::vector<int>& recvInd = commMeta.childRecvInd[proc];
            const int numRecvInd = recvInd.size();
            for( int k=0; k<numRecvInd; ++k )
            {
                const int iFrontLoc = recvInd[k];
                const T* recvRow = &recvValues[k*width];
                T* WRow = W.Buffer( iFrontLoc, 0 );
                const int WLDim = W.LDim();
                for( int jFront=0; jFront<width; ++jFront )
                    WRow[jFront*WLDim] += recvRow[jFront];
            }
        }
        std::vector<T>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );

        // Store this node's portion of the result
        X.distNodes[s-1] = WT;
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work1d.Empty();
}

template<typename T> 
inline void DistLowerMultiplyTranspose
( Orientation orientation, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, 
  DistNodalMultiVec<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistLowerMultiplyTranspose");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = X.Width();
    if( L.frontType != SYMM_1D && L.frontType != LDL_1D )
        throw std::logic_error("This multiply mode is not yet implemented");

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    const SymmFront<T>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        Matrix<T>& XRoot = X.localNodes.back();
        localRootFront.work = XRoot;
        FrontLowerMultiply
        ( orientation, diagOffset, localRootFront.frontL, XRoot );
    }
    else
    {
        const DistSymmFront<T>& rootFront = L.distFronts.back();
        DistMatrix<T,VC,STAR>& XRoot = X.distNodes.back();
        rootFront.work1d = XRoot;
        FrontLowerMultiply
        ( orientation, diagOffset, rootFront.front1dL, XRoot );
    }

    for( int s=numDistNodes-2; s>=0; --s )
    {
        const DistSymmNodeInfo& parentNode = info.distNodes[s+1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<T>& parentFront = L.distFronts[s+1];
        const DistSymmFront<T>& front = L.distFronts[s];
        const Grid& grid = front.front1dL.Grid();
        const Grid& parentGrid = parentFront.front1dL.Grid();
        mpi::Comm comm = grid.VCComm(); 
        mpi::Comm parentComm = parentGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int parentCommSize = mpi::CommSize( parentComm );

        // Set up a copy of the RHS in our workspace.
        DistMatrix<T,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( front.front1dL.Height(), width );
        DistMatrix<T,VC,STAR> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, node.size );
        Matrix<T>& localXT = 
          ( s>0 ? X.distNodes[s-1].Matrix() : X.localNodes.back() );
        WT.Matrix() = localXT;

        //
        // Set the bottom from the parent's workspace
        //

        // Pack the relevant portions of the parent's RHS's
        // (which are stored in 'work1d')
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
        std::vector<T> sendBuffer( sendBufferSize );

        DistMatrix<T,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            T* sendValues = &sendBuffer[sendDispls[proc]];
            const std::vector<int>& recvInd = commMeta.childRecvInd[proc];
            const int numRecvInd = recvInd.size();
            for( int k=0; k<numRecvInd; ++k )
            {
                const int iFrontLoc = recvInd[k];
                T* packedRow = &sendValues[k*width];
                const T* workRow = parentWork.LockedBuffer( iFrontLoc, 0 );
                const int workLDim = parentWork.LDim();
                for( int jFront=0; jFront<width; ++jFront )
                    packedRow[jFront] = workRow[jFront*workLDim];
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
        std::vector<T> recvBuffer( recvBufferSize );
#ifndef RELEASE
        VerifySendsAndRecvs( sendCounts, recvCounts, parentComm );
#endif

        // AllToAll to send and recv parent updates
        SparseAllToAll
        ( sendBuffer, sendCounts, sendDispls,
          recvBuffer, recvCounts, recvDispls, parentComm );
        std::vector<T>().swap( sendBuffer );
        std::vector<int>().swap( sendCounts );
        std::vector<int>().swap( sendDispls );

        // Unpack the updates using the send approach from the forward solve
        const std::vector<int>& myRelInd = 
            ( node.onLeft ? parentNode.leftRelInd : parentNode.rightRelInd );
        const int colShift = WB.ColShift();
        const int localHeight = WB.LocalHeight();
        for( int iUpdateLoc=0; iUpdateLoc<localHeight; ++iUpdateLoc )
        {
            const int iUpdate = colShift + iUpdateLoc*commSize;
            const int startRank = myRelInd[iUpdate] % parentCommSize;
            const T* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int jUpdate=0; jUpdate<width; ++jUpdate )
                WB.SetLocal(iUpdateLoc,jUpdate,recvBuf[jUpdate]);
            recvDispls[startRank] += width;
        }
        std::vector<T>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );

        // Make a copy of the unmodified RHS
        DistMatrix<T,VC,STAR> XNode( W );

        // Perform the multiply for this front
        if( s > 0 )
            FrontLowerMultiply
            ( orientation, diagOffset, front.front1dL, XNode );
        else
        {
            localRootFront.work = W.Matrix();
            FrontLowerMultiply
            ( orientation, diagOffset, 
              localRootFront.frontL, XNode.Matrix() );
        }

        // Store the supernode portion of the result
        DistMatrix<T,VC,STAR> XNodeT(grid), XNodeB(grid);
        PartitionDown
        ( XNode, XNodeT,
                 XNodeB, node.size );
        localXT = XNodeT.Matrix();
        XNode.Empty();
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERMULTIPLY_DIST_HPP
