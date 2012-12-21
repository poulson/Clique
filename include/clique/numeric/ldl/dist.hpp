/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F> 
void DistLDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistLDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL )
{
#ifndef RELEASE
    PushCallStack("DistLDL");
#endif
    // The bottom front is already computed, so just view it
    LocalSymmFront<F>& topLocalFront = L.localFronts.back();
    DistSymmFront<F>& bottomDistFront = L.distFronts[0];
    const Grid& bottomGrid = *info.distNodes[0].grid;
    bottomDistFront.front2dL.Empty();
    bottomDistFront.front2dL.LockedView
    ( topLocalFront.frontL.Height(), topLocalFront.frontL.Width(), 0, 0, 
      topLocalFront.frontL.LockedBuffer(), topLocalFront.frontL.LDim(), 
      bottomGrid );
    bottomDistFront.work2d.Empty();
    bottomDistFront.work2d.LockedView
    ( topLocalFront.work.Height(), topLocalFront.work.Width(), 0, 0,
      topLocalFront.work.LockedBuffer(), topLocalFront.work.LDim(),
      bottomGrid );

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

        const bool computeFactRecvIndices = 
            ( node.childFactRecvIndices.size() == 0 );

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
        const DistMatrix<F>& childUpdate = childFront.work2d;
        const bool isLeftChild = childNode.onLeft;
        std::vector<int> sendCounts(commSize), sendDispls(commSize);
        int sendBufferSize = 0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const int sendSize = node.numChildFactSendIndices[proc];
            sendCounts[proc] = sendSize;
            sendDispls[proc] = sendBufferSize;
            sendBufferSize += sendSize;
        }
        std::vector<F> sendBuffer( sendBufferSize );

        const std::vector<int>& myChildRelIndices = 
            ( isLeftChild ? node.leftChildRelIndices
                          : node.rightChildRelIndices );
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
                    childUpdate.GetLocal(iChildLocal,jChildLocal);
            }
        }
#ifndef RELEASE
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            if( packOffsets[proc]-sendDispls[proc] != 
                node.numChildFactSendIndices[proc] )
                throw std::logic_error("Error in packing stage");
        }
#endif
        packOffsets.clear();
        childFront.work2d.Empty();
        if( s == 1 )
            topLocalFront.work.Empty();

        // Set up the recv buffer for the AllToAll
        if( computeFactRecvIndices )
            ComputeFactRecvIndices( node, childNode );
        std::vector<int> recvCounts(commSize), recvDispls(commSize);
        int recvBufferSize=0;
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const int recvSize = node.childFactRecvIndices[proc].size()/2;
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

        // Unpack the child udpates (with an Axpy)
        front.work2d.SetGrid( front.front2dL.Grid() );
        front.work2d.Align( node.size % gridHeight, node.size % gridWidth );
        elem::Zeros( updateSize, updateSize, front.work2d );
        const int leftLocalWidth = front.front2dL.LocalWidth();
        const int topLocalHeight = 
            LocalLength<int>( node.size, grid.MCRank(), gridHeight );
        for( unsigned proc=0; proc<commSize; ++proc )
        {
            const F* recvValues = &recvBuffer[recvDispls[proc]];
            const std::deque<int>& recvIndices = 
                node.childFactRecvIndices[proc];
            const int numRecvIndexPairs = recvIndices.size()/2;
            for( int k=0; k<numRecvIndexPairs; ++k )
            {
                const int iFrontLocal = recvIndices[2*k+0];
                const int jFrontLocal = recvIndices[2*k+1];
                const F value = recvValues[k];
                if( jFrontLocal < leftLocalWidth )
                    front.front2dL.UpdateLocal
                    ( iFrontLocal, jFrontLocal, value );
                else
                    front.work2d.UpdateLocal
                    ( iFrontLocal-topLocalHeight, 
                      jFrontLocal-leftLocalWidth, value );
            }
        }
        recvBuffer.clear();
        recvCounts.clear();
        recvDispls.clear();
        if( computeFactRecvIndices )
            node.childFactRecvIndices.clear();

        // Now that the frontal matrix is set up, perform the factorization
        if( !blockLDL )
        {
            if( L.isHermitian )
                DistFrontLDL( ADJOINT, front.front2dL, front.work2d );
            else
                DistFrontLDL( TRANSPOSE, front.front2dL, front.work2d );

            // Store the diagonal in a [VC,* ] distribution
            DistMatrix<F,MD,STAR> diag( grid );
            front.front2dL.GetDiagonal( diag );
            front.diag.SetGrid( grid );
            front.diag = diag;
        }
        else
        {
            if( L.isHermitian )
                DistFrontBlockLDL( ADJOINT, front.front2dL, front.work2d );
            else
                DistFrontBlockLDL( TRANSPOSE, front.front2dL, front.work2d );
        }
    }
    L.localFronts.back().work.Empty();
    L.distFronts.back().work2d.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
