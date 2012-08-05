/*
   Clique: a scalable implementation of the multifrontal algorithm

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

namespace cliq {

template<typename F> 
void DistLowerForwardSolve
( UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

template<typename F>
void DistLowerBackwardSolve
( Orientation orientation, UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistLowerForwardSolve
( UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("DistLowerForwardSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = localX.Width();
    const SymmFrontType frontType = L.frontType;
    const bool frontsAre1d = FrontsAre1d( frontType );
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D );
    if( frontType != LDL_1D && 
        frontType != LDL_SELINV_1D && 
        frontType != LDL_SELINV_2D && 
        frontType != BLOCK_LDL_2D )
        throw std::logic_error("This solve mode is not yet implemented");
#ifndef RELEASE
    if( blockLDL && diag == UNIT )
        throw std::logic_error("Unit diagonal is nonsensical for block LDL");
#endif

    // Copy the information from the local portion into the distributed leaf
    const LocalSymmFront<F>& localRootFront = L.localFronts.back();
    const DistSymmFront<F>& distLeafFront = L.distFronts[0];
    const Grid& leafGrid = ( frontsAre1d ? distLeafFront.front1dL.Grid() 
                                         : distLeafFront.front2dL.Grid() );
    distLeafFront.work1d.LockedView
    ( localRootFront.work.Height(), localRootFront.work.Width(), 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim(), 
      leafGrid );
    
    // Perform the distributed portion of the forward solve
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& childFront = L.distFronts[s-1];
        const DistSymmFront<F>& front = L.distFronts[s];
        const Grid& childGrid = ( frontsAre1d ? childFront.front1dL.Grid()
                                              : childFront.front2dL.Grid() );
        const Grid& grid = ( frontsAre1d ? front.front1dL.Grid()
                                         : front.front2dL.Grid() );
        mpi::Comm comm = grid.VCComm();
        mpi::Comm childComm = childGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int childCommSize = mpi::CommSize( childComm );
        const int frontHeight = ( frontsAre1d ? front.front1dL.Height()
                                              : front.front2dL.Height() );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( frontHeight, width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.View( localX, node.localOffset1d, 0, node.localSize1d, width );
        WT.LocalMatrix() = localXT;
        elem::MakeZeros( WB );

        // Pack our child's update
        DistMatrix<F,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childNode.size;
        DistMatrix<F,VC,STAR> childUpdate;
        childUpdate.LockedView( childW, childNode.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int sendSize = node.numChildSolveSendIndices[proc]*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const bool isLeftChild = childNode.onLeft;
        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? node.leftChildRelIndices
                          : node.rightChildRelIndices );
        const int updateColShift = childUpdate.ColShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        std::vector<int> packOffsets = sendDispls;
        for( int iChildLocal=0; iChildLocal<updateLocalHeight; ++iChildLocal )
        {
            const int iChild = updateColShift + iChildLocal*childCommSize;
            const int destRank = myChildRelIndices[iChild] % commSize;
            F* packBuf = &sendBuffer[packOffsets[destRank]];
            for( int jChild=0; jChild<width; ++jChild )
                packBuf[jChild] = childUpdate.GetLocal(iChildLocal,jChild);
            packOffsets[destRank] += width;
        }
        packOffsets.clear();
        childW.Empty();
        if( s == 1 )
            L.localFronts.back().work.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int recvSize = node.childSolveRecvIndices[proc].size()*width;
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
        sendBuffer.clear();
        sendCounts.clear();
        sendDispls.clear();

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = 
                node.childSolveRecvIndices[proc];
            for( unsigned k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                const F* recvRow = &recvValues[k*width];
                F* WRow = W.LocalBuffer( iFrontLocal, 0 );
                const int WLDim = W.LocalLDim();
                for( int j=0; j<width; ++j )
                    WRow[j*WLDim] += recvRow[j];
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Now that the RHS is set up, perform this node's solve
        if( frontType == LDL_1D )
            DistFrontLowerForwardSolve( diag, front.front1dL, W );
        else if( frontType == LDL_SELINV_1D )
            DistFrontFastLowerForwardSolve( diag, front.front1dL, W );
        else if( frontType == LDL_SELINV_2D )
            DistFrontFastLowerForwardSolve( diag, front.front2dL, W );
        else // frontType == BLOCK_LDL_2D
            DistFrontBlockLowerForwardSolve( front.front2dL, W );

        // Store this node's portion of the result
        localXT = WT.LocalMatrix();
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work1d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void DistLowerBackwardSolve
( Orientation orientation, UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("DistLowerBackwardSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = localX.Width();
    const SymmFrontType frontType = L.frontType;
    const bool frontsAre1d = FrontsAre1d( frontType );
    const bool blockLDL = ( frontType == BLOCK_LDL_2D );
    if( frontType != LDL_1D && 
        frontType != LDL_SELINV_1D && 
        frontType != LDL_SELINV_2D && 
        frontType != BLOCK_LDL_2D )
        throw std::logic_error("This solve mode is not yet implemented");
#ifndef RELEASE
    if( blockLDL && diag == UNIT )
        throw std::logic_error("Unit diagonal is nonsensical for block LDL");
#endif

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    const LocalSymmFront<F>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        localRootFront.work.View
        ( rootNode.size, width, 
          localX.Buffer(rootNode.localOffset1d,0), localX.LDim() );
        if( !blockLDL )
            LocalFrontLowerBackwardSolve
            ( orientation, diag, localRootFront.frontL, localRootFront.work );
        else
            LocalFrontBlockLowerBackwardSolve
            ( orientation, localRootFront.frontL, localRootFront.work );
    }
    else
    {
        const DistSymmFront<F>& rootFront = L.distFronts.back();
        const Grid& rootGrid = ( frontsAre1d ? rootFront.front1dL.Grid() 
                                             : rootFront.front2dL.Grid() );
        rootFront.work1d.View
        ( rootNode.size, width, 0,
          localX.Buffer(rootNode.localOffset1d,0), localX.LDim(), rootGrid );
        if( frontType == LDL_1D )
            DistFrontLowerBackwardSolve
            ( orientation, diag, rootFront.front1dL, rootFront.work1d );
        else if( frontType == LDL_SELINV_1D )
            DistFrontFastLowerBackwardSolve
            ( orientation, diag, rootFront.front1dL, rootFront.work1d );
        else if( frontType == LDL_SELINV_2D )
            DistFrontFastLowerBackwardSolve
            ( orientation, diag, rootFront.front2dL, rootFront.work1d );
        else
            DistFrontBlockLowerBackwardSolve
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
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.View( localX, node.localOffset1d, 0, node.localSize1d, width );
        WT.LocalMatrix() = localXT;

        //
        // Set the bottom from the parent
        //

        // Pack the updates using the recv approach from the forward solve
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int sendSize = 
                parentNode.childSolveRecvIndices[proc].size()*width;
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendValues = &sendBuffer[sendDispls[proc]];
            const std::deque<int>& recvIndices = 
                parentNode.childSolveRecvIndices[proc];
            for( unsigned k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                F* sendRow = &sendValues[k*width];
                const F* workRow = 
                    parentWork.LockedLocalBuffer( iFrontLocal, 0 );
                const int workLDim = parentWork.LocalLDim();
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
            const int recvSize = 
                parentNode.numChildSolveSendIndices[proc]*width;
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
        sendBuffer.clear();
        sendCounts.clear();
        sendDispls.clear();

        // Unpack the updates using the send approach from the forward solve
        const bool isLeftChild = node.onLeft;
        const std::vector<int>& myRelIndices = 
            ( isLeftChild ? parentNode.leftChildRelIndices
                          : parentNode.rightChildRelIndices );
        const int updateColShift = WB.ColShift();
        const int updateLocalHeight = WB.LocalHeight();
        for( int iUpdateLocal=0; 
                 iUpdateLocal<updateLocalHeight; ++iUpdateLocal )
        {
            const int iUpdate = updateColShift + iUpdateLocal*commSize;
            const int startRank = myRelIndices[iUpdate] % parentCommSize;
            const F* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int j=0; j<width; ++j )
                WB.SetLocal(iUpdateLocal,j,recvBuf[j]);
            recvDispls[startRank] += width;
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Call the custom node backward solve
        if( s > 0 )
        {
            if( frontType == LDL_1D )
                DistFrontLowerBackwardSolve
                ( orientation, diag, front.front1dL, W );
            else if( frontType == LDL_SELINV_1D )
                DistFrontFastLowerBackwardSolve
                ( orientation, diag, front.front1dL, W );
            else if( frontType == LDL_SELINV_2D )
                DistFrontFastLowerBackwardSolve
                ( orientation, diag, front.front2dL, W );
            else // frontType == BLOCK_LDL_2D
                DistFrontBlockLowerBackwardSolve
                ( orientation, front.front2dL, front.work1d );
        }
        else
        {
            localRootFront.work.View( W.LocalMatrix() );
            if( !blockLDL )
                LocalFrontLowerBackwardSolve
                ( orientation, diag, localRootFront.frontL, 
                  localRootFront.work );
            else
                LocalFrontBlockLowerBackwardSolve
                ( orientation, localRootFront.frontL, localRootFront.work );
        }

        // Store this node's portion of the result
        localXT = WT.LocalMatrix();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
