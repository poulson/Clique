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
DistNodalMultiVec<F>::DistNodalMultiVec()
: height_(0), width_(0), localHeight_(0)
{ }

template<typename F>
inline
DistNodalMultiVec<F>::DistNodalMultiVec
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMultiVec::DistNodalMultiVec");
#endif
    Pull( inverseMap, info, X );
}

template<typename F>
inline void
DistNodalMultiVec<F>::Pull
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMultiVec::Pull");
#endif
    height_ = X.Height();
    width_ = X.Width();

    // Traverse our part of the elimination tree to see how many indices we need
    int numRecvIndices=0;
    const int numLocal = info.localNodes.size();
    for( int s=0; s<numLocal; ++s )
        numRecvIndices += info.localNodes[s].size;
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
        numRecvIndices += info.distNodes[s].solveMeta1d.localSize;
    localHeight_ = numRecvIndices;
    
    // Fill the set of indices that we need to map to the original ordering
    int offset=0;
    std::vector<int> mappedIndices( numRecvIndices );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
            mappedIndices[offset++] = node.offset+t;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();
        const int alignment = 0;
        const int shift = Shift( gridRank, alignment, gridSize );
        for( int t=shift; t<node.size; t+=gridSize )
            mappedIndices[offset++] = node.offset+t;
    }
#ifndef RELEASE
    if( offset != numRecvIndices )
        throw std::logic_error("mappedIndices was filled incorrectly");
#endif

    // Convert the indices to the original ordering
    inverseMap.Translate( mappedIndices );

    // Figure out how many entries each process owns that we need
    mpi::Comm comm = X.Comm();
    const int commSize = mpi::CommSize( comm );
    std::vector<int> recvSizes( commSize, 0 );
    const int blocksize = X.Blocksize();
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++recvSizes[q];
    }
    std::vector<int> recvOffsets( commSize );
    offset=0;
    for( int q=0; q<commSize; ++q )
    {
        recvOffsets[q] = offset;
        offset += recvSizes[q];
    }
    std::vector<int> recvIndices( numRecvIndices );
    std::vector<int> offsets = recvOffsets;
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        recvIndices[offsets[q]++] = i;
    }

    // Coordinate for the coming AllToAll to exchange the indices of X
    std::vector<int> sendSizes( commSize );
    mpi::AllToAll( &recvSizes[0], 1, &sendSizes[0], 1, comm );
    int numSendIndices=0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSendIndices;
        numSendIndices += sendSizes[q];
    }

    // Request the indices
    std::vector<int> sendIndices( numSendIndices );
    mpi::AllToAll
    ( &recvIndices[0], &recvSizes[0], &recvOffsets[0],
      &sendIndices[0], &sendSizes[0], &sendOffsets[0], comm );

    // Fulfill the requests
    std::vector<F> sendValues( numSendIndices*width_ );
    const int firstLocalRow = X.FirstLocalRow();
    for( int s=0; s<numSendIndices; ++s )
        for( int j=0; j<width_; ++j )
            sendValues[s*width_+j] = 
                X.GetLocal( sendIndices[s]-firstLocalRow, j );

    // Reply with the values
    std::vector<F> recvValues( numRecvIndices*width_ );
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] *= width_;
        sendOffsets[q] *= width_;
        recvSizes[q] *= width_;
        recvOffsets[q] *= width_;
    }
    mpi::AllToAll
    ( &sendValues[0], &sendSizes[0], &sendOffsets[0],
      &recvValues[0], &recvSizes[0], &recvOffsets[0], comm );
    std::vector<F>().swap( sendValues );
    std::vector<int>().swap( sendSizes );
    std::vector<int>().swap( sendOffsets );

    // Unpack the values
    offset = 0;
    offsets = recvOffsets;
    localNodes.resize( numLocal );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        localNodes[s].ResizeTo( node.size, width_ );
        for( int t=0; t<node.size; ++t )
        {
            const int i = mappedIndices[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width_; ++j )
                localNodes[s].Set( t, j, recvValues[offsets[q]++] );
        }
    }
    distNodes.resize( numDist-1 );
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();
        const int alignment = 0;
        const int shift = Shift( gridRank, alignment, gridSize );
        const int localHeight = Length( node.size, shift, gridSize );
        distNodes[s-1].ResizeTo( localHeight, width_ );
        for( int tLoc=0; tLoc<localHeight; ++tLoc )
        {
            const int i = mappedIndices[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width_; ++j )
                distNodes[s-1].Set( tLoc, j, recvValues[offsets[q]++] );
        }
    }
#ifndef RELEASE
    if( offset != numRecvIndices )
        throw std::logic_error("Unpacked wrong number of indices");
#endif
}

template<typename F>
inline void
DistNodalMultiVec<F>::Push
( const DistMap& inverseMap, const DistSymmInfo& info,
        DistMultiVec<F>& X ) const
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMultiVec::Push");
#endif
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    mpi::Comm comm = rootNode.comm;
    const int height = rootNode.size + rootNode.offset;
    const int width = Width();
    X.SetComm( comm );
    X.ResizeTo( height, width );

    const int commSize = mpi::CommSize( comm );
    const int blocksize = X.Blocksize();
    const int localHeight = X.LocalHeight();
    const int firstLocalRow = X.FirstLocalRow();
    const int numDist = info.distNodes.size();
    const int numLocal = info.localNodes.size();

    // Fill the set of indices that we need to map to the original ordering
    const int numSendIndices = LocalHeight();
    int offset=0;
    std::vector<int> mappedIndices( numSendIndices );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
            mappedIndices[offset++] = node.offset+t;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();
        const int alignment = 0;
        const int shift = Shift( gridRank, alignment, gridSize );
        for( int t=shift; t<node.size; t+=gridSize )
            mappedIndices[offset++] = node.offset+t;
    }

    // Convert the indices to the original ordering
    inverseMap.Translate( mappedIndices );

    // Figure out how many indices each process owns that we need to send
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numSendIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++sendSizes[q];
    }
    std::vector<int> sendOffsets( commSize );
    offset=0;
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = offset;
        offset += sendSizes[q];
    }

    // Pack the send indices and values
    offset=0;
    std::vector<F> sendValues( numSendIndices*width );
    std::vector<int> sendIndices( numSendIndices );
    std::vector<int> offsets = sendOffsets;
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
        {
            const int i = mappedIndices[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width; ++j )
                sendValues[offsets[q]*width+j] = localNodes[s].Get(t,j);    
            sendIndices[offsets[q]++] = i;
        }
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& grid = *node.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();
        const int alignment = 0;
        const int shift = Shift( gridRank, alignment, gridSize );
        const int localHeight = Length( node.size, shift, gridSize );

        for( int tLoc=0; tLoc<localHeight; ++tLoc )
        {
            const int i = mappedIndices[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width; ++j )
                sendValues[offsets[q]*width+j] = distNodes[s-1].Get(tLoc,j);
            sendIndices[offsets[q]++] = i;
        }
    }

    // Coordinate for the coming AllToAll to exchange the indices of x
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm );
    int numRecvIndices=0;
    std::vector<int> recvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvOffsets[q] = numRecvIndices;
        numRecvIndices += recvSizes[q];
    }
#ifndef RELEASE
    if( numRecvIndices != localHeight )
        throw std::logic_error("numRecvIndices was not equal to local height");
#endif

    // Send the indices
    std::vector<int> recvIndices( numRecvIndices );
    mpi::AllToAll
    ( &sendIndices[0], &sendSizes[0], &sendOffsets[0],
      &recvIndices[0], &recvSizes[0], &recvOffsets[0], comm );

    // Send the values
    std::vector<F> recvValues( numRecvIndices*width );
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] *= width;
        sendOffsets[q] *= width;
        recvSizes[q] *= width;
        recvOffsets[q] *= width;
    }
    mpi::AllToAll
    ( &sendValues[0], &sendSizes[0], &sendOffsets[0],
      &recvValues[0], &recvSizes[0], &recvOffsets[0], comm );
    std::vector<F>().swap( sendValues );
    std::vector<int>().swap( sendSizes );
    std::vector<int>().swap( sendOffsets );

    // Unpack the values
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = recvIndices[s];
        const int iLocal = i - firstLocalRow;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= localHeight )
            throw std::logic_error("iLocal was out of bounds");
#endif
        for( int j=0; j<width; ++j )
            X.SetLocal( iLocal, j, recvValues[s*width+j] );
    }
}

template<typename F>
inline int
DistNodalMultiVec<F>::Height() const
{ return height_; }

template<typename F>
inline int
DistNodalMultiVec<F>::Width() const
{ return width_; }

template<typename F>
inline int
DistNodalMultiVec<F>::LocalHeight() const
{ return localHeight_; }

} // namespace cliq
