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
    if( numSupernodes == 0 )
        return;

    // The bottom front is already computed, so just view it
    distL.fronts[0].LocalMatrix().LockedView( localL.fronts.back() );

    // Perform the distributed portion of the factorization
    std::vector<int>::const_iterator it;
    for( unsigned k=1; k<numSupernodes; ++k )
    {
        const symbolic::DistSymmFactSupernode& childSymbSN = S.supernodes[k-1];
        const symbolic::DistSymmFactSupernode& symbSN = S.supernodes[k];
        const DistMatrix<F,MC,MR>& childFront = distL.fronts[k-1];
        DistMatrix<F,MC,MR>& front = distL.fronts[k];
        const bool computeRecvIndices = ( symbSN.childRecvIndices.size() == 0 );

        // Grab this front's grid information
        const Grid& grid = front.Grid();
        mpi::Comm comm = grid.VCComm();
        const unsigned commRank = mpi::CommRank( comm );
        const unsigned commSize = mpi::CommSize( comm );
        const unsigned gridHeight = grid.Height();
        const unsigned gridWidth = grid.Width();

        // Grab the child's grid information
        const Grid& childGrid = childFront.Grid();
        mpi::Comm childComm = childGrid.VCComm();
        const unsigned childCommRank = mpi::CommRank( childComm );
        const unsigned childCommSize = mpi::CommSize( childComm );
        const unsigned childGridHeight = childGrid.Height();
        const unsigned childGridWidth = childGrid.Width();
        const unsigned childGridRow = childGrid.MCRank();
        const unsigned childGridCol = childGrid.MRRank();

#ifndef RELEASE
        if( front.Height() != symbSN.size+symbSN.lowerStruct.size() ||
            front.Width()  != symbSN.size+symbSN.lowerStruct.size() )
            throw std::logic_error("Front was not the proper size");
#endif

        // Pack our child's updates
        DistMatrix<F,MC,MR> childUpdate(childGrid);
        const int updateSize = childFront.Height()-childSymbSN.size;
        childUpdate.LockedView
        ( childFront, 
          childSymbSN.size, childSymbSN.size, updateSize, updateSize );
        const bool isLeftChild = ( commRank < commSize/2 );
        it = std::max_element
             ( symbSN.numChildSendIndices.begin(), 
               symbSN.numChildSendIndices.end() );
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
        if( computeRecvIndices )
            symbolic::ComputeRecvIndices( symbSN );
        int recvPortionSize = mpi::MIN_COLL_MSG;
        for( int i=0; i<commSize; ++i )
        {
            const int thisPortion = symbSN.childRecvIndices[i].size();
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
            const std::deque<int>& recvIndices = symbSN.childRecvIndices[proc];
            for( int k=0; k<recvIndices.size(); ++k )
            {
                const int iFrontLocal = recvIndices[2*k+0];
                const int jFrontLocal = recvIndices[2*k+1];
                const F value = recvValues[k];
                front.UpdateLocalEntry( iFrontLocal, jFrontLocal, value );
            }
        }
        recvBuffer.clear();
        if( computeRecvIndices )
            symbSN.childRecvIndices.clear();

        // Now that the frontal matrix is set up, perform the factorization
        DistSupernodeLDL( orientation, front, symbSN.size );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<float>& localL,
        numeric::DistSymmFact<float>& distL );

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<double>& localL,
        numeric::DistSymmFact<double>& distL );

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
        numeric::DistSymmFact<std::complex<float> >& distL );

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
        numeric::DistSymmFact<std::complex<double> >& distL );

