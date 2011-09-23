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

// This routine could be modified later so that it uses much less memory
// by replacing the '=' redistributions with piece-by-piece redistributions.
template<typename F>
void clique::numeric::SetSolveMode( DistSymmFact<F>& distL, SolveMode mode )
{
#ifndef RELEASE
    PushCallStack("numeric::SetSolveMode");
#endif
    // Check if this call can be a no-op
    if( mode == distL.mode ) 
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    distL.mode = mode;
    const int numSupernodes = distL.supernodes.size();    
    if( numSupernodes == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    DistSymmFactSupernode<F>& leafSN = distL.supernodes[0];
    if( mode == FEW_RHS )
    {
        leafSN.front1d.LocalMatrix().View( leafSN.front2d.LocalMatrix() );
        for( int k=1; k<numSupernodes; ++k )
        {
            DistSymmFactSupernode<F>& sn = distL.supernodes[k];
            sn.front1d = sn.front2d;
            sn.front2d.Empty();
        }
    }
    else
    {
        leafSN.front2d.LocalMatrix().View( leafSN.front1d.LocalMatrix() );
        for( int k=1; k<numSupernodes; ++k )
        {
            DistSymmFactSupernode<F>& sn = distL.supernodes[k];
            sn.front2d = sn.front1d;
            sn.front1d.Empty();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S, // can't be const due to map...
  const numeric::LocalSymmFact<F>& localL,
        numeric::DistSymmFact<F>& distL )
{
#ifndef RELEASE
    PushCallStack("numeric::DistLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = S.supernodes.size();
    distL.mode = MANY_RHS;
    if( numSupernodes == 0 )
        return;

    // The bottom front is already computed, so just view it
    const LocalSymmFactSupernode<F>& topLocalSN = localL.supernodes.back();
    DistSymmFactSupernode<F>& bottomDistSN = distL.supernodes[0];
    const Grid& bottomGrid = bottomDistSN.front2d.Grid();
    bottomDistSN.front2d.Empty(); // eventually this can be removed...
    bottomDistSN.front2d.LockedView
    ( topLocalSN.front.Height(), topLocalSN.front.Width(), 0, 0, 
      topLocalSN.front.LockedBuffer(), topLocalSN.front.LDim(), bottomGrid );

    // Perform the distributed portion of the factorization
    std::vector<int>::const_iterator it;
    for( unsigned k=1; k<numSupernodes; ++k )
    {
        const symbolic::DistSymmFactSupernode& childSymbSN = S.supernodes[k-1];
        const symbolic::DistSymmFactSupernode& symbSN = S.supernodes[k];
        const DistSymmFactSupernode<F>& childNumSN = distL.supernodes[k-1];
        DistSymmFactSupernode<F>& numSN = distL.supernodes[k];

        const bool computeFactRecvIndices = 
            ( symbSN.childFactRecvIndices.size() == 0 );

        // Grab this front's grid information
        const Grid& grid = numSN.front2d.Grid();
        mpi::Comm comm = grid.VCComm();
        const unsigned commRank = mpi::CommRank( comm );
        const unsigned commSize = mpi::CommSize( comm );
        const unsigned gridHeight = grid.Height();
        const unsigned gridWidth = grid.Width();

        // Grab the child's grid information
        const Grid& childGrid = childNumSN.front2d.Grid();
        mpi::Comm childComm = childGrid.VCComm();
        const unsigned childCommRank = mpi::CommRank( childComm );
        const unsigned childCommSize = mpi::CommSize( childComm );
        const unsigned childGridHeight = childGrid.Height();
        const unsigned childGridWidth = childGrid.Width();
        const unsigned childGridRow = childGrid.MCRank();
        const unsigned childGridCol = childGrid.MRRank();

#ifndef RELEASE
        if( numSN.front2d.Height() != symbSN.size+symbSN.lowerStruct.size() ||
            numSN.front2d.Width()  != symbSN.size+symbSN.lowerStruct.size() )
            throw std::logic_error("Front was not the proper size");
#endif

        // Pack our child's update
        DistMatrix<F,MC,MR> childUpdate(childGrid);
        const int updateSize = childNumSN.front2d.Height()-childSymbSN.size;
        childUpdate.LockedView
        ( childNumSN.front2d, 
          childSymbSN.size, childSymbSN.size, updateSize, updateSize );
        const bool isLeftChild = ( commRank < commSize/2 );
        it = std::max_element
             ( symbSN.numChildFactSendIndices.begin(), 
               symbSN.numChildFactSendIndices.end() );
        const int sendPortionSize = std::max(*it,mpi::MIN_COLL_MSG);
        std::vector<F> sendBuffer( sendPortionSize*commSize );

        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? symbSN.leftChildRelIndices
                          : symbSN.rightChildRelIndices );
        const int updateRowAlignment = childUpdate.RowAlignment();
        const int updateColShift = childUpdate.ColShift();
        const int updateRowShift = childUpdate.RowShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        const int updateLocalWidth = childUpdate.LocalWidth();
        // Initialize the offsets to each process's chunk
        std::vector<int> sendOffsets( commSize );
        for( int proc=0; proc<commSize; ++proc )
            sendOffsets[proc] = proc*sendPortionSize;
        for( int jChildLocal=0; jChildLocal<updateLocalWidth; ++jChildLocal )
        {
            const int jChild = updateRowShift + jChildLocal*childGridWidth;
            const int destGridCol = myChildRelIndices[jChild] % gridWidth;

            const int align = (jChild+updateRowAlignment) % childGridHeight;
            const int shift = 
                (childGridRow+childGridHeight-align) % childGridHeight;
            const int localColShift = 
                (jChild+shift-updateColShift) / childGridHeight;
            for( int iChildLocal=localColShift; 
                     iChildLocal<updateLocalHeight; ++iChildLocal )
            {
                const int iChild = updateColShift + iChildLocal*childGridHeight;
                const int destGridRow = myChildRelIndices[iChild] % gridHeight;

                const int destRank = destGridRow + destGridCol*gridHeight;
                sendBuffer[sendOffsets[destRank]++] = 
                    childUpdate.GetLocalEntry(iChildLocal,jChildLocal);
            }
        }
        // Reset the offsets to their original values
        for( int proc=0; proc<commSize; ++proc )
            sendOffsets[proc] = proc*sendPortionSize;

        // AllToAll to send and receive the child updates
        if( computeFactRecvIndices )
            symbolic::ComputeFactRecvIndices( symbSN );
        int recvPortionSize = mpi::MIN_COLL_MSG;
        for( int proc=0; proc<commSize; ++proc )
        {
            const int thisPortion = symbSN.childFactRecvIndices[proc].size();
            recvPortionSize = std::max(thisPortion,recvPortionSize);
        }
        std::vector<F> recvBuffer( recvPortionSize*commSize );
        mpi::AllToAll
        ( &sendBuffer[0], sendPortionSize, 
          &recvBuffer[0], recvPortionSize, comm );
        sendBuffer.clear();

        // Unpack the child udpates (with an Axpy)
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[proc*recvPortionSize];
            const std::deque<int>& recvIndices = 
                symbSN.childFactRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[2*k+0];
                const int jFrontLocal = recvIndices[2*k+1];
                const F value = recvValues[k];
                numSN.front2d.UpdateLocalEntry
                ( iFrontLocal, jFrontLocal, value );
            }
        }
        recvBuffer.clear();
        if( computeFactRecvIndices )
            symbSN.childFactRecvIndices.clear();

        // Now that the frontal matrix is set up, perform the factorization
        DistSupernodeLDL( orientation, numSN.front2d, symbSN.size );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::SetSolveMode
( DistSymmFact<float>& distL, SolveMode mode );
template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<float>& localL,
        numeric::DistSymmFact<float>& distL );

template void clique::numeric::SetSolveMode
( DistSymmFact<double>& distL, SolveMode mode );
template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<double>& localL,
        numeric::DistSymmFact<double>& distL );

template void clique::numeric::SetSolveMode
( DistSymmFact<std::complex<float> >& distL, SolveMode mode );
template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
        numeric::DistSymmFact<std::complex<float> >& distL );

template void clique::numeric::SetSolveMode
( DistSymmFact<std::complex<double> >& distL, SolveMode mode );
template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
        numeric::DistSymmFact<std::complex<double> >& distL );
