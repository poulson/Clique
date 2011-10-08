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
#include "clique.hpp"
using namespace elemental;

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L, 
        Matrix<F>& localX )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistLDLForwardSolve");
#endif
    const int numSupernodes = S.dist.supernodes.size();
    const int width = localX.Width();
    if( L.dist.mode == MANY_RHS )
        throw std::logic_error("This solve mode is not yet implemented");
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Copy the information from the local portion into the distributed leaf
    const LocalSymmFront<F>& localRootFront = L.local.fronts.back();
    const DistSymmFront<F>& distLeafFront = L.dist.fronts[0];
    distLeafFront.work1d.LockedView
    ( localRootFront.work.Height(), localRootFront.work.Width(), 0,
      localRootFront.work.LockedBuffer(), localRootFront.work.LDim(), 
      distLeafFront.front1d.Grid() );
    
    // Perform the distributed portion of the forward solve
    for( int s=1; s<numSupernodes; ++s )
    {
        const DistSymmFactSupernode& childSN = S.dist.supernodes[s-1];
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const DistSymmFront<F>& childFront = L.dist.fronts[s-1];
        const DistSymmFront<F>& front = L.dist.fronts[s];
        const Grid& childGrid = childFront.front1d.Grid();
        const Grid& grid = front.front1d.Grid();
        mpi::Comm comm = grid.VCComm();
        mpi::Comm childComm = childGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int commRank = mpi::CommRank( comm );
        const int childCommSize = mpi::CommSize( childComm );
        const int childCommRank = mpi::CommRank( childComm );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( front.front1d.Height(), width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.View( localX, sn.localOffset1d, 0, sn.localSize1d, width );
        WT.LocalMatrix() = localXT;
        WB.SetToZero();

        // Pack our child's update
        DistMatrix<F,VC,STAR>& childW = childFront.work1d;
        const int updateSize = childW.Height()-childSN.size;
        DistMatrix<F,VC,STAR> childUpdate;
        childUpdate.LockedView( childW, childSN.size, 0, updateSize, width );
        int sendBufferSize = 0;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int actualSend = sn.numChildSolveSendIndices[proc]*width;
            const int thisSend = std::max(actualSend,1);
            sendCounts[proc] = thisSend;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += thisSend;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const bool isLeftChild = ( commRank < commSize/2 );
        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? sn.leftChildRelIndices
                          : sn.rightChildRelIndices );
        const int updateColShift = childUpdate.ColShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        std::vector<int> packOffsets = sendDispls;
        for( int iChildLocal=0; iChildLocal<updateLocalHeight; ++iChildLocal )
        {
            const int iChild = updateColShift + iChildLocal*childCommSize;
            const int destRank = myChildRelIndices[iChild] % commSize;
            F* packBuf = &sendBuffer[packOffsets[destRank]];
            for( int jChild=0; jChild<width; ++jChild )
                packBuf[jChild] = childUpdate.GetLocalEntry(iChildLocal,jChild);
            packOffsets[destRank] += width;
        }
        packOffsets.clear();
        childW.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        for( int proc=0; proc<commSize; ++proc )
        {
            const int actualRecv = sn.childSolveRecvIndices[proc].size()*width;
            const int thisRecv = std::max(actualRecv,1);
            recvCounts[proc] = thisRecv;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += thisRecv;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        // Verify the send and recv counts match
        std::vector<int> actualRecvCounts(commSize);
        mpi::AllToAll
        ( &sendCounts[0],       1,
          &actualRecvCounts[0], 1, comm );
        for( int proc=0; proc<commSize; ++proc )
        {
            if( actualRecvCounts[proc] != recvCounts[proc] )
            {
                std::ostringstream msg;
                msg << "Expected recv count of " << recvCounts[proc]
                    << " but recv'd " << actualRecvCounts[proc]
                    << " from process " << proc << " for supernode "
                    << s << "\n";
                throw std::logic_error( msg.str().c_str() );
            }
        }
        actualRecvCounts.clear();
#endif

        // AllToAll to send and receive the child updates
        mpi::AllToAll
        ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
          &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm );
        sendBuffer.clear();
        sendCounts.clear();
        sendDispls.clear();

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = sn.childSolveRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                const F* recvRow = &recvValues[k*width];
                F* WRow = W.LocalBuffer( iFrontLocal, 0 );
                const int WLDim = W.LocalLDim();
                for( int jFront=0; jFront<width; ++jFront )
                    WRow[jFront*WLDim] += recvRow[jFront];
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Now that the RHS is set up, perform this supernode's solve
        DistFrontLDLForwardSolve( sn.size, front.front1d, W );

        // Store the supernode portion of the result
        localXT = WT.LocalMatrix();
    }

    // Free the distributed and local root work buffers
    L.local.fronts.back().work.Empty();
    L.dist.fronts.back().work1d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F representa a real or complex field
void clique::numeric::DistLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX, bool checkIfSingular )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistLDLDiagonalSolve");
#endif
    const int numSupernodes = S.dist.supernodes.size();
    const int width = localX.Width();

    for( int s=1; s<numSupernodes; ++s )
    {
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const DistSymmFront<F>& front = L.dist.fronts[s];

        Matrix<F> localXT;
        localXT.View( localX, sn.localOffset1d, 0, sn.localSize1d, width );

        // Solve against the s'th supernode using the front
        DistMatrix<F,VC,STAR> FTL;
        FTL.LockedView( front.front1d, 0, 0, sn.size, sn.size );
        DistMatrix<F,VC,STAR> dTL;
        FTL.GetDiagonal( dTL );
        basic::DiagonalSolve
        ( LEFT, NORMAL, dTL.LockedLocalMatrix(), localXT, checkIfSingular );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& localX )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistLDLBackwardSolve");
#endif
    const int numSupernodes = S.dist.supernodes.size();
    const int width = localX.Width();
    if( L.dist.mode == MANY_RHS )
        throw std::logic_error("This solve mode is not yet implemented");
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Directly operate on the root separator's portion of the right-hand sides
    const DistSymmFactSupernode& rootSN = S.dist.supernodes.back();
    const DistSymmFront<F>& rootFront = L.dist.fronts.back();
    const Grid& rootGrid = rootFront.front1d.Grid();
    rootFront.work1d.View
    ( rootSN.size, width, 0,
      localX.Buffer(rootSN.localOffset1d,0), localX.LDim(), rootGrid );
    DistFrontLDLBackwardSolve
    ( orientation, rootSN.size, rootFront.front1d, rootFront.work1d );

    std::vector<int>::const_iterator it;
    for( int s=numSupernodes-2; s>=0; --s )
    {
        const DistSymmFactSupernode& parentSN = S.dist.supernodes[s+1];
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        const DistSymmFront<F>& parentFront = L.dist.fronts[s+1];
        const DistSymmFront<F>& front = L.dist.fronts[s];
        const Grid& grid = front.front1d.Grid();
        const Grid& parentGrid = parentFront.front1d.Grid();
        mpi::Comm comm = grid.VCComm(); 
        mpi::Comm parentComm = parentGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int commRank = mpi::CommRank( comm );
        const int parentCommSize = mpi::CommSize( parentComm );
        const int parentCommRank = mpi::CommRank( parentComm );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = front.work1d;
        W.SetGrid( grid );
        W.ResizeTo( front.front1d.Height(), width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.View( localX, sn.localOffset1d, 0, sn.localSize1d, width );
        WT.LocalMatrix() = localXT;

        //
        // Set the bottom from the parent
        //

        // Pack the updates using the recv approach from the forward solve
        int sendBufferSize = 0;
        std::vector<int> sendCounts(parentCommSize), sendDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int actualSend = 
                parentSN.childSolveRecvIndices[proc].size()*width;
            const int thisSend = std::max(actualSend,1);
            sendCounts[proc] = thisSend;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += thisSend;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        DistMatrix<F,VC,STAR>& parentWork = parentFront.work1d;
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            F* sendValues = &sendBuffer[sendDispls[proc]];
            const std::deque<int>& recvIndices = 
                parentSN.childSolveRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                F* sendRow = &sendValues[k*width];
                const F* workRow = parentWork.LocalBuffer( iFrontLocal, 0 );
                const int workLDim = parentWork.LocalLDim();
                for( int jFront=0; jFront<width; ++jFront )
                    sendRow[jFront] = workRow[jFront*workLDim];
            }
        }
        parentWork.Empty();

        // Set up the receive buffer
        int recvBufferSize = 0;
        std::vector<int> recvCounts(parentCommSize), recvDispls(parentCommSize);
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            const int actualRecv = 
                parentSN.numChildSolveSendIndices[proc]*width;
            const int thisRecv = std::max(actualRecv,1);
            recvCounts[proc] = thisRecv;
            recvDispls[proc] = recvBufferSize;
            recvBufferSize += thisRecv;
        }
        std::vector<F> recvBuffer( recvBufferSize );
#ifndef RELEASE
        // Verify the send and recv counts match
        std::vector<int> actualRecvCounts(parentCommSize);
        mpi::AllToAll
        ( &sendCounts[0],       1,
          &actualRecvCounts[0], 1, parentComm );
        for( int proc=0; proc<parentCommSize; ++proc )
        {
            if( actualRecvCounts[proc] != recvCounts[proc] )
            {
                std::ostringstream msg;
                msg << "Expected recv count of " << recvCounts[proc]
                    << " but recv'd " << actualRecvCounts[proc]
                    << " from process " << proc << " for supernode "
                    << s << "\n";
                throw std::logic_error( msg.str().c_str() );
            }
        }
        actualRecvCounts.clear();
#endif

        // AllToAll to send and recv parent updates
        mpi::AllToAll
        ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
          &recvBuffer[0], &recvCounts[0], &recvDispls[0], parentComm );
        sendBuffer.clear();
        sendCounts.clear();
        sendDispls.clear();

        // Unpack the updates using the send approach from the forward solve
        const bool isLeftChild = ( parentCommRank < parentCommSize/2 );
        const std::vector<int>& myRelIndices = 
            ( isLeftChild ? parentSN.leftChildRelIndices
                          : parentSN.rightChildRelIndices );
        const int updateColShift = WB.ColShift();
        const int updateLocalHeight = WB.LocalHeight();
        for( int iUpdateLocal=0; 
                 iUpdateLocal<updateLocalHeight; ++iUpdateLocal )
        {
            const int iUpdate = updateColShift + iUpdateLocal*commSize;
            const int startRank = myRelIndices[iUpdate] % parentCommSize;
            const F* recvBuf = &recvBuffer[recvDispls[startRank]];
            for( int jUpdate=0; jUpdate<width; ++jUpdate )
                WB.SetLocalEntry(iUpdateLocal,jUpdate,recvBuf[jUpdate]);
            recvDispls[startRank] += width;
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();

        // Call the custom supernode backward solve
        DistFrontLDLBackwardSolve( orientation, sn.size, front.front1d, W );

        // Store the supernode portion of the result
        localXT = WT.LocalMatrix();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& localX );
template void clique::numeric::DistLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& localX,
        bool checkIfSingular );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& localX );
template void clique::numeric::DistLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& localX,
        bool checkIfSingular );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& localX );
template void clique::numeric::DistLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& localX,
        bool checkIfSingular );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& localX );
template void clique::numeric::DistLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& localX,
        bool checkIfSingular );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& localX );
