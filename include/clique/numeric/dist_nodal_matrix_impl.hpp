/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix()
: height_(0), width_(0)
{ }

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::DistNodalMatrix");
#endif
    Pull( inverseMap, info, X );
}

template<typename F>
inline void
DistNodalMatrix<F>::Pull
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::Pull");
#endif
    mpi::Comm comm = X.Comm();
    const int commSize = mpi::CommSize( comm );
    const int height = X.Height();
    const int width = X.Width();
    const int blocksize = X.Blocksize();
    const int firstLocalRow = X.FirstLocalRow();
    const int numDist = info.distNodes.size();
    const int numLocal = info.localNodes.size();

    height_ = height;
    width_ = width;

    // Traverse our part of the elimination tree to see which row indices
    // we will need
    int numRecvRowInd=0;
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        numRecvRowInd += node.size;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        numRecvRowInd += node.solveMeta2d.localHeight;
    }

    // Fill the set of indices that we need to map to the original ordering
    int offset=0;
    std::vector<int> mappedInd( numRecvRowInd );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
            mappedInd[offset++] = node.offset+t;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridRow = grid.Row();
        const int gridHeight = grid.Height();
        const int colAlign = 0;
        const int colShift = Shift( gridRow, colAlign, gridHeight );
        for( int t=colShift; t<node.size; t+=gridHeight )
            mappedInd[offset++] = node.offset+t;
    }
#ifndef RELEASE
    if( offset != numRecvRowInd )
        throw std::logic_error("mappedInd was filled incorrectly");
#endif

    // Convert the row indices to the original ordering
    inverseMap.Translate( mappedInd );

    // Count the number of indices that we will request from each process
    offset = 0;
    std::vector<int> numRecvInd( commSize, 0 );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            numRecvInd[q] += width;
        }
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridRow = grid.Row();
        const int gridCol = grid.Col();
        const int gridHeight = grid.Height();
        const int gridWidth = grid.Width();
        const int colAlign = 0;
        const int rowAlign = 0;
        const int colShift = Shift( gridRow, colAlign, gridHeight );
        const int localWidth = Length( width, gridCol, rowAlign, gridWidth );
        for( int t=colShift; t<node.size; t+=gridHeight )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            numRecvInd[q] += localWidth;
        }
    }

    // Fill the set of indices that we will request from each process
    int totalRecvInd=0;
    std::vector<int> recvSizes( commSize ), recvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        totalRecvInd += numRecvInd[q];
        recvSizes[q] = 2*numRecvInd[q];
        recvOffsets[q] = 2*totalRecvInd;
    }
    offset = 0;
    std::vector<int> recvInd( 2*totalRecvInd ); 
    std::vector<int> offsets = recvOffsets;
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int u=0; u<width; ++u )
            {
                recvInd[offsets[q]++] = i;      
                recvInd[offsets[q]++] = u;
            }
        }
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridRow = grid.Row();
        const int gridCol = grid.Col();
        const int gridHeight = grid.Height();
        const int gridWidth = grid.Width();
        const int colAlign = 0;
        const int rowAlign = 0;
        const int colShift = Shift( gridRow, colAlign, gridHeight );
        const int rowShift = Shift( gridCol, rowAlign, gridWidth );
        for( int t=colShift; t<node.size; t+=gridHeight )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int u=rowShift; u<width; u+=gridWidth )
            {
                recvInd[offsets[q]++] = i;
                recvInd[offsets[q]++] = u;
            }
        }
    }

    // Coordinate for the coming AllToAll to exchange the indices of X
    std::vector<int> numSendInd( commSize );
    mpi::AllToAll( &numRecvInd[0], 1, &numSendInd[0], 1, comm );
    int totalSendInd=0;
    std::vector<int> sendSizes( commSize ), sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        totalSendInd += numSendInd[q];
        sendSizes[q] = 2*numSendInd[q];
        sendOffsets[q] = 2*totalSendInd;
    }

    // Request the indices
    std::vector<int> sendInd( 2*totalSendInd );
    mpi::AllToAll
    ( &recvInd[0], &recvSizes[0], &recvOffsets[0],
      &sendInd[0], &sendSizes[0], &sendOffsets[0], comm );
    // TODO? swap out recvInd?

    // Fulfill the requests
    std::vector<F> sendValues( totalSendInd );
    for( int s=0; s<totalSendInd; ++s )
        sendValues[s] = 
            X.GetLocal( sendInd[2*s]-firstLocalRow, sendInd[2*s+1] );

    // Reply with the values
    std::vector<F> recvValues( totalRecvInd );
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] /= 2;
        sendOffsets[q] /= 2;
        recvSizes[q] /= 2;
        recvOffsets[q] /= 2;
    }
    mpi::AllToAll
    ( &sendValues[0], &sendSizes[0], &sendOffsets[0],
      &recvValues[0], &recvSizes[0], &recvOffsets[0], comm );
    std::vector<F>().swap( sendValues );
    std::vector<int>().swap( sendSizes );
    std::vector<int>().swap( sendOffsets );

    // Unpack the values
    // HERE: Dang...need to reorganize data structure
    throw std::logic_error("This routine is not yet finished");
    /*
    offset=0;
    offsets = recvOffsets;
    matrix.ResizeTo( numRecvInd, width );
    */
}

template<typename F>
inline void
DistNodalMatrix<F>::Push
( const DistMap& inverseMap, const DistSymmInfo& info,
        DistMultiVec<F>& X ) const
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::Push");
#endif
    throw std::logic_error("This must be rewritten");
}

template<typename F>
inline int
DistNodalMatrix<F>::Height() const
{ return height_; }

template<typename F>
inline int
DistNodalMatrix<F>::Width() const
{ return width_; }

} // namespace cliq
