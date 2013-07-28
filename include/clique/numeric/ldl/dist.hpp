/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LDL_DIST_HPP
#define CLIQ_NUMERIC_LDL_DIST_HPP

namespace cliq {

template<typename F> 
void 
DistLDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void 
DistLDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL )
{
#ifndef RELEASE
    CallStackEntry entry("DistLDL");
#endif
    // The bottom front is already computed, so just view it
    SymmFront<F>& topLocalFront = L.localFronts.back();
    DistSymmFront<F>& bottomDistFront = L.distFronts[0];
    const Grid& bottomGrid = *info.distNodes[0].grid;
    bottomDistFront.front2dL.Empty();
    bottomDistFront.front2dL.LockedAttach
    ( topLocalFront.frontL.Height(), topLocalFront.frontL.Width(), 0, 0, 
      topLocalFront.frontL.LockedBuffer(), topLocalFront.frontL.LDim(), 
      bottomGrid );
    bottomDistFront.work2d.Empty();
    bottomDistFront.work2d.LockedAttach
    ( topLocalFront.work.Height(), topLocalFront.work.Width(), 0, 0,
      topLocalFront.work.LockedBuffer(), topLocalFront.work.LDim(),
      bottomGrid );
    if( !blockLDL )
    {
        // Store the diagonal in a [VC,* ] distribution
        DistMatrix<F,MD,STAR> diag( bottomGrid );
        bottomDistFront.front2dL.GetDiagonal( diag );
        bottomDistFront.diag1d.SetGrid( bottomGrid );
        bottomDistFront.diag1d = diag;
        elem::SetDiagonal( bottomDistFront.front2dL, F(1) );
    }

    // Perform the distributed portion of the factorization
    const unsigned numDistNodes = info.distNodes.size();
    for( unsigned s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const int updateSize = node.lowerStruct.size();
        DistSymmFront<F>& childFront = L.distFronts[s-1];
        DistSymmFront<F>& front = L.distFronts[s];
        front.work2d.Empty();
#ifndef RELEASE
        if( front.front2dL.Height() != node.size+updateSize ||
            front.front2dL.Width() != node.size )
            throw std::logic_error("Front was not the proper size");
#endif

        // Grab this front's grid information
        const Grid& grid = front.front2dL.Grid();
        mpi::Comm comm = grid.VCComm();
        const unsigned commSize = mpi::CommSize( comm );
        const unsigned gridHeight = grid.Height();
        const unsigned gridWidth = grid.Width();

        // Grab the child's grid information
        const Grid& childGrid = childFront.front2dL.Grid();
        const unsigned childGridHeight = childGrid.Height();
        const unsigned childGridWidth = childGrid.Width();

        // Pack our child's update
        const FactorCommMeta& commMeta = node.factorMeta;
        const DistMatrix<F>& childUpdate = childFront.work2d;
        const bool onLeft = childNode.onLeft;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        int sendBufferSize = 0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const int sendSize = commMeta.numChildSendInds[proc];
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const std::vector<int>& myChildRelInds = 
            ( onLeft ? node.leftRelInds : node.rightRelInds );
        const int updateColShift = childUpdate.ColShift();
        const int updateRowShift = childUpdate.RowShift();
        const int updateLocalHeight = childUpdate.LocalHeight();
        const int updateLocalWidth = childUpdate.LocalWidth();
        std::vector<int> packOffs = sendDispls;
        for( int jChildLoc=0; jChildLoc<updateLocalWidth; ++jChildLoc )
        {
            const int jChild = updateRowShift + jChildLoc*childGridWidth;
            const int destGridCol = myChildRelInds[jChild] % gridWidth;
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
                const int iChild = updateColShift + iChildLoc*childGridHeight;
                if( iChild >= jChild )
                {
                    const int destGridRow = myChildRelInds[iChild] % gridHeight;
                    const int destRank = destGridRow + destGridCol*gridHeight;
                    sendBuffer[packOffs[destRank]++] = 
                        childUpdate.GetLocal(iChildLoc,jChildLoc);
                }
            }
        }
#ifndef RELEASE
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            if( packOffs[proc]-sendDispls[proc] != 
                commMeta.numChildSendInds[proc] )
                throw std::logic_error("Error in packing stage");
        }
#endif
        std::vector<int>().swap( packOffs );
        childFront.work2d.Empty();
        if( s == 1 )
            topLocalFront.work.Empty();

        // Set up the recv buffer for the AllToAll
        const bool computeFactRecvInds = ( commMeta.childRecvInds.size() == 0 );
        if( computeFactRecvInds )
            ComputeFactRecvInds( node, childNode );
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        int recvBufferSize=0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const int recvSize = commMeta.childRecvInds[proc].size()/2;
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

        // Unpack the child udpates (with an Axpy)
        front.work2d.SetGrid( front.front2dL.Grid() );
        front.work2d.Align( node.size % gridHeight, node.size % gridWidth );
        Zeros( front.work2d, updateSize, updateSize );
        const int leftLocalWidth = front.front2dL.LocalWidth();
        const int topLocalHeight = 
            Length<int>( node.size, grid.Row(), gridHeight );
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const F* recvVals = &recvBuffer[recvDispls[proc]];
            const std::vector<int>& recvInds = commMeta.childRecvInds[proc];
            const int numRecvIndPairs = recvInds.size()/2;
            for( int k=0; k<numRecvIndPairs; ++k )
            {
                const int iFrontLoc = recvInds[2*k+0];
                const int jFrontLoc = recvInds[2*k+1];
                const F value = recvVals[k];
#ifndef RELEASE
                const int iFront = grid.Row() + iFrontLoc*gridHeight;
                const int jFront = grid.Col() + jFrontLoc*gridWidth;
                if( iFront < jFront )
                    throw std::logic_error("Tried to update upper triangle");
#endif
                if( jFrontLoc < leftLocalWidth )
                    front.front2dL.UpdateLocal( iFrontLoc, jFrontLoc, value );
                else
                    front.work2d.UpdateLocal
                    ( iFrontLoc-topLocalHeight, jFrontLoc-leftLocalWidth, 
                      value );
            }
        }
        std::vector<F>().swap( recvBuffer );
        std::vector<int>().swap( recvCounts );
        std::vector<int>().swap( recvDispls );
        if( computeFactRecvInds )
            commMeta.EmptyChildRecvIndices();

        // Now that the frontal matrix is set up, perform the factorization
        if( !blockLDL )
        {
            if( L.isHermitian )
                FrontLDL( ADJOINT, front.front2dL, front.work2d );
            else
                FrontLDL( TRANSPOSE, front.front2dL, front.work2d );

            // Store the diagonal in a [VC,* ] distribution
            DistMatrix<F,MD,STAR> diag( grid );
            front.front2dL.GetDiagonal( diag );
            front.diag1d.SetGrid( grid );
            front.diag1d = diag;
            elem::SetDiagonal( front.front2dL, F(1) );
        }
        else
        {
            if( L.isHermitian )
                FrontBlockLDL( ADJOINT, front.front2dL, front.work2d );
            else
                FrontBlockLDL( TRANSPOSE, front.front2dL, front.work2d );
        }
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work2d.Empty();
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_DIST_HPP
