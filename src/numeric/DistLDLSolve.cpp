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
( const symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<F>& localL,
  const numeric::DistSymmFact<F>& distL,
        Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("numeric::DistLDLForwardSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = localX.Width();
    if( distL.mode == MANY_RHS )
        throw std::logic_error("This solve mode is not yet implemented");
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Copy the information from the local portion into the distributed root
    const LocalSymmFactSupernode<F>& topLocalSN = localL.supernodes.back();
    const DistSymmFactSupernode<F>& bottomDistSN = distL.supernodes[0];
    bottomDistSN.work2d.LocalMatrix().LockedView( topLocalSN.work );
    
    // Perform the distributed portion of the forward solve
    std::vector<int>::const_iterator it;
    for( int k=1; k<numSupernodes; ++k )
    {
        const symbolic::DistSymmFactSupernode& childSymbSN = S.supernodes[k-1];
        const symbolic::DistSymmFactSupernode& symbSN = S.supernodes[k];
        const DistSymmFactSupernode<F>& childNumSN = distL.supernodes[k-1];
        const DistSymmFactSupernode<F>& numSN = distL.supernodes[k];
        const Grid& childGrid = childNumSN.front1d.Grid();
        const Grid& grid = numSN.front1d.Grid();
        mpi::Comm comm = grid.VCComm();
        mpi::Comm childComm = childGrid.VCComm();
        const int commSize = mpi::CommSize( comm );
        const int commRank = mpi::CommRank( comm );
        const int childCommSize = mpi::CommSize( childComm );
        const int childCommRank = mpi::CommRank( childComm );

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = numSN.work1d;
        W.ResizeTo( numSN.front1d.Height(), width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        WT.View( W, 0, 0, symbSN.size, width );
        WB.View( W, symbSN.size, 0, W.Height()-symbSN.size, width );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.LockedView
        ( localX, symbSN.localOffset1d, 0, symbSN.localSize1d, width );
        WT.LocalMatrix() = localXT;
        WB.SetToZero();

        // Pack our child's update
        DistMatrix<F,VC,STAR> childUpdate(childGrid);
        const int updateSize = childNumSN.work1d.Height()-childSymbSN.size;
        childUpdate.LockedView
        ( childNumSN.work1d, childSymbSN.size, 0, updateSize, width );
        const bool isLeftChild = ( commRank < commSize/2 );
        it = std::max_element
             ( symbSN.numChildSolveSendIndices.begin(), 
               symbSN.numChildSolveSendIndices.end() );
        const int sendPortionSize = std::max((*it)*width,mpi::MIN_COLL_MSG);
        std::vector<F> sendBuffer( sendPortionSize*commSize );

        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? symbSN.leftChildRelIndices
                          : symbSN.rightChildRelIndices );
        const int updateColShift = childUpdate.ColShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        // Initialize the offsets to each process's chunk
        std::vector<int> sendOffsets( commSize );
        for( int proc=0; proc<commSize; ++proc )
            sendOffsets[proc] = proc*sendPortionSize;
        for( int iChildLocal=0; iChildLocal<updateLocalHeight; ++iChildLocal )
        {
            const int iChild = updateColShift + iChildLocal*childCommSize;
            const int destRank = myChildRelIndices[iChild] % commSize;
            F* sendBuf = &sendBuffer[sendOffsets[destRank]];
            for( int jChild=0; jChild<width; ++jChild )
                sendBuf[jChild] = childUpdate.GetLocalEntry(iChildLocal,jChild);
            sendOffsets[destRank] += width;
        }
        // Free the child work buffer
        childNumSN.work1d.Empty();
        // Reset the offsets to their original values
        for( int proc=0; proc<commSize; ++proc )
            sendOffsets[proc] = proc*sendPortionSize;

        // AllToAll to send and receive the child updates
        int recvPortionSize = mpi::MIN_COLL_MSG;
        for( int proc=0; proc<commSize; ++proc )
        {
            const int thisColumn = symbSN.childSolveRecvIndices[proc].size();
            recvPortionSize = std::max(thisColumn*width,recvPortionSize); 
        }
        std::vector<F> recvBuffer( recvPortionSize*commSize );
        mpi::AllToAll
        ( &sendBuffer[0], sendPortionSize,
          &recvBuffer[0], recvPortionSize, comm );
        sendBuffer.clear();

        // Unpack the child updates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[proc*recvPortionSize];
            const std::deque<int>& recvIndices = 
                symbSN.childSolveRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[k];
                const F* recvRow = &recvValues[k*width];
                F* workRow = numSN.work1d.LocalBuffer( iFrontLocal, 0 );
                const int workLDim = numSN.work1d.LocalLDim();
                for( int jFront=0; jFront<width; ++jFront )
                    workRow[jFront*workLDim] = recvRow[jFront];
            }
        }
        recvBuffer.clear();

        // Now that the RHS is set up, perform this supernode's solve
        DistSupernodeLDLForwardSolve( symbSN.size, numSN.front1d, W );

        // Store the supernode portion of the result
        localXT = WT.LocalMatrix();
    }

    // Free the distributed and local root work buffers
    localL.supernodes.back().work.Empty();
    distL.supernodes.back().work1d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::DistSymmFact& S,
  const numeric::DistSymmFact<F>& L,
        Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("numeric::DistLDLBackwardSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = localX.Width();
    if( L.mode == MANY_RHS )
        throw std::logic_error("This solve mode is not yet implemented");
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Directly operate on the root separator's portion of the right-hand sides
    const symbolic::DistSymmFactSupernode& rootSymbSN = S.supernodes.back();
    const DistSymmFactSupernode<F>& rootNumSN = L.supernodes.back();
    const Grid& rootGrid = rootNumSN.front1d.Grid();
    rootNumSN.work1d.View
    ( rootSymbSN.size, width, 0,
      localX.Buffer(rootSymbSN.localOffset1d,0), localX.LDim(), rootGrid );
    DistSupernodeLDLBackwardSolve
    ( orientation, rootSymbSN.size, rootNumSN.front1d, rootNumSN.work1d );

    for( int k=numSupernodes-2; k>=0; --k )
    {
        const symbolic::DistSymmFactSupernode& symbSN = S.supernodes[k];
        const numeric::DistSymmFactSupernode<F>& numSN = L.supernodes[k];

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = numSN.work1d;
        W.SetGrid( numSN.front1d.Grid() );
        W.ResizeTo( numSN.front1d.Height(), width );
        DistMatrix<F,VC,STAR> WT, WB;
        PartitionDown
        ( W, WT,
             WB, symbSN.size );

        // Pull in the relevant information from the RHS
        Matrix<F> localXT;
        localXT.LockedView
        ( localX, symbSN.localOffset1d, 0, symbSN.localSize1d, width );
        WT.LocalMatrix() = localXT;

        // Set the bottom from the parent
        DistMatrix<F,VC,STAR>& parentWork = L.supernodes[k+1].work1d;
        // HERE: 
        // Pack the updates using the recv approach from the forward solve
        // AllToAll
        // Unpack the updates using the send approach from the forward solve

        // Free the parent's workspace
        parentWork.Empty();

        // Call the custom supernode backward solve
        DistSupernodeLDLBackwardSolve
        ( orientation, symbSN.size, numSN.front1d, W );

        // Store the supernode portion of the result
        localXT = WT.LocalMatrix();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<float>& localL,
  const numeric::DistSymmFact<float>& distL,
        Matrix<float>& localX );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::DistSymmFact& S,
  const numeric::DistSymmFact<float>& L,
        Matrix<float>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<double>& localL,
  const numeric::DistSymmFact<double>& distL,
        Matrix<double>& localX );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::DistSymmFact& S,
  const numeric::DistSymmFact<double>& L,
        Matrix<double>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
  const numeric::DistSymmFact<std::complex<float> >& distL,
        Matrix<std::complex<float> >& localX );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::DistSymmFact& S,
  const numeric::DistSymmFact<std::complex<float> >& L,
        Matrix<std::complex<float> >& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
  const numeric::DistSymmFact<std::complex<double> >& distL,
        Matrix<std::complex<double> >& localX );
template void clique::numeric::DistLDLBackwardSolve
( Orientation orientation,
  const symbolic::DistSymmFact& S,
  const numeric::DistSymmFact<std::complex<double> >& L,
        Matrix<std::complex<double> >& localX );
