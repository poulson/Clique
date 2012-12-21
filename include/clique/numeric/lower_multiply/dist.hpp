/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T> 
void DistLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX );

template<typename T> 
void DistLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T> 
inline void DistLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )
{
#ifndef RELEASE
    PushCallStack("DistLowerMultiplyNormal");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = localX.Width();
    if( L.frontType != SYMM_1D && L.frontType != LDL_1D )
        throw std::logic_error("This multiply mode is not yet implemented");

    // Copy the information from the local portion into the distributed leaf
    const LocalSymmFront<T>& localRootFront = L.localFronts.back();
    const DistSymmFront<T>& distLeafFront = L.distFronts[0];
    distLeafFront.work1d.LockedView
    ( localRootFront.work.Height(), localRootFront.work.Width(), 0,
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
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<T> localXT;
        localXT.View( localX, node.localOffset1d, 0, node.localSize1d, width );
        WT.LocalMatrix() = localXT;
        elem::MakeZeros( WB );

        // Now that the right-hand side is set up, perform the multiply
        DistFrontLowerMultiply( NORMAL, diag, diagOffset, front.front1dL, W );

        // Pack our child's update
        DistMatrix<T,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childNode.size;
        DistMatrix<T,VC,STAR> childUpdate;
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
        std::vector<T> sendBuffer( sendBufferSize );

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
            T* packBuf = &sendBuffer[packOffsets[destRank]];
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
        std::vector<T> recvBuffer( recvBufferSize );
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
            const T* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = 
                node.childSolveRecvIndices[proc];
            const int numRecvIndices = recvIndices.size();
            for( int k=0; k<numRecvIndices; ++k )
            {
                const int iFrontLocal = recvIndices[k];
                const T* recvRow = &recvValues[k*width];
                T* WRow = W.LocalBuffer( iFrontLocal, 0 );
                const int WLDim = W.LocalLDim();
                for( int jFront=0; jFront<width; ++jFront )
                    WRow[jFront*WLDim] += recvRow[jFront];
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Store this node's portion of the result
        localXT = WT.LocalMatrix();
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work1d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void DistLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )
{
#ifndef RELEASE
    PushCallStack("DistLowerMultiplyTranspose");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = localX.Width();
    if( L.frontType != SYMM_1D && L.frontType != LDL_1D )
        throw std::logic_error("This multiply mode is not yet implemented");

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    const LocalSymmFront<T>& localRootFront = L.localFronts.back();
    if( numDistNodes == 1 )
    {
        Matrix<T> XRoot;
        XRoot.View
        ( rootNode.size, width, 
          localX.Buffer(rootNode.localOffset1d,0), localX.LDim() );
        localRootFront.work = XRoot;
        LocalFrontLowerMultiply
        ( orientation, diag, diagOffset, localRootFront.frontL, XRoot );
    }
    else
    {
        const DistSymmFront<T>& rootFront = L.distFronts.back();
        const Grid& rootGrid = rootFront.front1dL.Grid();
        DistMatrix<T,VC,STAR> XRoot(rootGrid);
        XRoot.View
        ( rootNode.size, width, 0,
          localX.Buffer(rootNode.localOffset1d,0), localX.LDim(), rootGrid );
        rootFront.work1d = XRoot; // store the RHS for use by the children
        DistFrontLowerMultiply
        ( orientation, diag, diagOffset, rootFront.front1dL, XRoot );
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
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<T> localXT;
        localXT.View( localX, node.localOffset1d, 0, node.localSize1d, width );
        WT.LocalMatrix() = localXT;

        //
        // Set the bottom from the parent's workspace
        //

        // Pack the relevant portions of the parent's RHS's
        // (which are stored in 'work1d')
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
        std::vector<T> sendBuffer( sendBufferSize );

        DistMatrix<T,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            T* sendValues = &sendBuffer[sendDispls[proc]];
            const std::deque<int>& recvIndices = 
                parentNode.childSolveRecvIndices[proc];
            const int numRecvIndices = recvIndices.size();
            for( int k=0; k<numRecvIndices; ++k )
            {
                const int iFrontLocal = recvIndices[k];
                T* packedRow = &sendValues[k*width];
                const T* workRow = 
                    parentWork.LockedLocalBuffer( iFrontLocal, 0 );
                const int workLDim = parentWork.LocalLDim();
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
            const int recvSize = 
                parentNode.numChildSolveSendIndices[proc]*width;
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
            const T* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int jUpdate=0; jUpdate<width; ++jUpdate )
                WB.SetLocal(iUpdateLocal,jUpdate,recvBuf[jUpdate]);
            recvDispls[startRank] += width;
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Make a copy of the unmodified RHS
        DistMatrix<T,VC,STAR> XNode( W );

        // Perform the multiply for this front
        if( s > 0 )
            DistFrontLowerMultiply
            ( orientation, diag, diagOffset, front.front1dL, XNode );
        else
        {
            localRootFront.work = W.LocalMatrix();
            LocalFrontLowerMultiply
            ( orientation, diag, diagOffset, 
              localRootFront.frontL, XNode.LocalMatrix() );
        }

        // Store the supernode portion of the result
        DistMatrix<T,VC,STAR> XNodeT(grid), XNodeB(grid);
        elem::PartitionDown
        ( XNode, XNodeT,
                 XNodeB, node.size );
        localXT = XNodeT.LocalMatrix();
        XNode.Empty();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
