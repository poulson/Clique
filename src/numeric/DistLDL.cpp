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

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDL
( Orientation orientation, symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::DistLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numDistSupernodes = S.dist.supernodes.size();
    L.dist.mode = MANY_RHS;

    // The bottom front is already computed, so just view it
    LocalSymmFront<F>& topLocalFront = L.local.fronts.back();
    DistSymmFront<F>& bottomDistFront = L.dist.fronts[0];
    const Grid& bottomGrid = *S.dist.supernodes[0].grid;
    bottomDistFront.front2dL.Empty();
    bottomDistFront.front2dL.LockedView
    ( topLocalFront.frontL.Height(), topLocalFront.frontL.Width(), 0, 0, 
      topLocalFront.frontL.LockedBuffer(), topLocalFront.frontL.LDim(), 
      bottomGrid );
    bottomDistFront.front2dR.Empty();
    bottomDistFront.front2dR.LockedView
    ( topLocalFront.frontR.Height(), topLocalFront.frontR.Width(), 0, 0,
      topLocalFront.frontR.LockedBuffer(), topLocalFront.frontR.LDim(),
      bottomGrid );

    // Perform the distributed portion of the factorization
    std::vector<int>::const_iterator it;
    for( unsigned s=1; s<numDistSupernodes; ++s )
    {
        const DistSymmFactSupernode& childSN = S.dist.supernodes[s-1];
        const DistSymmFactSupernode& sn = S.dist.supernodes[s];
        DistSymmFront<F>& childFront = L.dist.fronts[s-1];
        DistSymmFront<F>& front = L.dist.fronts[s];

        const bool computeFactRecvIndices = 
            ( sn.childFactRecvIndices.size() == 0 );

        // Grab this front's grid information
        const Grid& grid = front.front2dL.Grid();
        mpi::Comm comm = grid.VCComm();
        const unsigned commRank = mpi::CommRank( comm );
        const unsigned commSize = mpi::CommSize( comm );
        const unsigned gridHeight = grid.Height();
        const unsigned gridWidth = grid.Width();

        // Grab the child's grid information
        const Grid& childGrid = childFront.front2dL.Grid();
        mpi::Comm childComm = childGrid.VCComm();
        const unsigned childCommRank = mpi::CommRank( childComm );
        const unsigned childCommSize = mpi::CommSize( childComm );
        const unsigned childGridHeight = childGrid.Height();
        const unsigned childGridWidth = childGrid.Width();
        const unsigned childGridRow = childGrid.MCRank();
        const unsigned childGridCol = childGrid.MRRank();

#ifndef RELEASE
        if( front.front2dL.Height() != sn.size+sn.lowerStruct.size() ||
            front.front2dL.Width() != sn.size )
            throw std::logic_error("Front was not the proper size");
#endif
        front.front2dR.SetGrid( front.front2dL.Grid() );
        front.front2dR.Align( 0, sn.size%grid.Width() );
        front.front2dR.ResizeTo
        ( front.front2dL.Height(), sn.lowerStruct.size() );
        front.front2dR.SetToZero();

        // Pack our child's update
        DistMatrix<F,MC,MR> childUpdate;
        const int updateSize = childFront.front2dR.Width();
        childUpdate.LockedView
        ( childFront.front2dR, 
          childSN.size, 0, updateSize, updateSize );
        const bool isLeftChild = ( commRank < commSize/2 );
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        int sendBufferSize = 0;
        for( int proc=0; proc<commSize; ++proc )
        {
            const int actualSend = sn.numChildFactSendIndices[proc];
            const int thisSend = std::max(actualSend,1);
            sendCounts[proc] = thisSend;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += thisSend;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? sn.leftChildRelIndices
                          : sn.rightChildRelIndices );
        const int updateRowAlignment = childUpdate.RowAlignment();
        const int updateColShift = childUpdate.ColShift();
        const int updateRowShift = childUpdate.RowShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        const int updateLocalWidth = childUpdate.LocalWidth();
        std::vector<int> packOffsets = sendDispls;
        for( int jChildLocal=0; jChildLocal<updateLocalWidth; ++jChildLocal )
        {
            const int jChild = updateRowShift + jChildLocal*childGridWidth;
            const int destGridCol = myChildRelIndices[jChild] % gridWidth;

            int localColShift;
            if( updateColShift > jChild )
                localColShift = 0;
            else if( (jChild-updateColShift) % childGridHeight == 0 )
                localColShift = (jChild-updateColShift)/childGridHeight;
            else
                localColShift = (jChild-updateColShift)/childGridHeight + 1;
            for( int iChildLocal=localColShift; 
                     iChildLocal<updateLocalHeight; ++iChildLocal )
            {
                const int iChild = updateColShift + iChildLocal*childGridHeight;
                const int destGridRow = myChildRelIndices[iChild] % gridHeight;

                const int destRank = destGridRow + destGridCol*gridHeight;
                sendBuffer[packOffsets[destRank]++] = 
                    childUpdate.GetLocalEntry(iChildLocal,jChildLocal);
            }
        }
#ifndef RELEASE
        for( int proc=0; proc<commSize; ++proc )
        {
            if( packOffsets[proc]-sendDispls[proc] != 
                sn.numChildFactSendIndices[proc] )
                throw std::logic_error("Error in packing stage");
        }
#endif
        packOffsets.clear();
        childFront.front2dR.Empty();
        if( s == 1 )
            topLocalFront.frontR.Empty();

        // Set up the recv buffer for the AllToAll
        if( computeFactRecvIndices )
            ComputeFactRecvIndices( sn, childSN );
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        int recvBufferSize = 0;
        for( int proc=0; proc<commSize; ++proc )
        {
            const int actualRecv = sn.childFactRecvIndices[proc].size()/2;
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

        // Unpack the child udpates (with an Axpy)
        const int leftLocalWidth = front.front2dL.LocalWidth();
        for( int proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = sn.childFactRecvIndices[proc];
            const int numRecvIndexPairs = recvIndices.size()/2;
            for( int k=0; k<numRecvIndexPairs; ++k )
            {
                const int iFrontLocal = recvIndices[2*k+0];
                const int jFrontLocal = recvIndices[2*k+1];
                const F value = recvValues[k];
                if( jFrontLocal < leftLocalWidth )
                    front.front2dL.UpdateLocalEntry
                    ( iFrontLocal, jFrontLocal, value );
                else
                    front.front2dR.UpdateLocalEntry
                    ( iFrontLocal, jFrontLocal-leftLocalWidth, value );
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();
        if( computeFactRecvIndices )
            sn.childFactRecvIndices.clear();

        // Now that the frontal matrix is set up, perform the factorization
        DistFrontLDL( orientation, front.front2dL, front.front2dR );
    }
    L.local.fronts.back().frontR.Empty();
    L.dist.fronts.back().front2dR.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<float>& L );

template void clique::numeric::DistLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<double>& L );

template void clique::numeric::DistLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<std::complex<float> >& L );

template void clique::numeric::DistLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<std::complex<double> >& L );
