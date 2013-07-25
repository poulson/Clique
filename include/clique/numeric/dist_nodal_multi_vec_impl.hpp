/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_DISTNODALMULTIVEC_IMPL_HPP
#define CLIQ_NUMERIC_DISTNODALMULTIVEC_IMPL_HPP

namespace cliq {

template<typename F>
inline
DistNodalMultiVec<F>::DistNodalMultiVec()
: height_(0), width_(0)
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
inline
DistNodalMultiVec<F>::DistNodalMultiVec( const DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMultiVec::DistNodalMultiVec");
#endif
    *this = X;
}

template<typename F>
inline const DistNodalMultiVec<F>&
DistNodalMultiVec<F>::operator=( const DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMultiVec::operator=");
#endif
    height_ = X.Height();
    width_ = X.Width();

    // Copy over the nontrivial distributed nodes
    const int numDist = X.distNodes.size();
    distNodes.resize( numDist );
    for( int s=0; s<numDist; ++s )
    {
        distNodes[s].SetGrid( X.distNodes[s].Grid() );
        distNodes[s] = X.distNodes[s];
    }

    // Copy over the local nodes
    const int numLocal = X.localNodes.size();
    localNodes.resize( numLocal );
    for( int s=0; s<numLocal; ++s )
        localNodes[s] = X.localNodes[s];

    return *this;
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
    int numRecvInd=0;
    const int numLocal = info.localNodes.size();
    for( int s=0; s<numLocal; ++s )
        numRecvInd += info.localNodes[s].size;
    const int numDist = info.distNodes.size();
    for( int s=1; s<numDist; ++s )
        numRecvInd += info.distNodes[s].multiVecMeta.localSize;
    
    // Fill the set of indices that we need to map to the original ordering
    int offset=0;
    std::vector<int> mappedInd( numRecvInd );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& nodeInfo = info.localNodes[s];
        for( int t=0; t<nodeInfo.size; ++t )
            mappedInd[offset++] = nodeInfo.offset+t;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        const Grid& grid = *nodeInfo.grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.VCRank();
        const int alignment = 0;
        const int shift = Shift( gridRank, alignment, gridSize );
        for( int t=shift; t<nodeInfo.size; t+=gridSize )
            mappedInd[offset++] = nodeInfo.offset+t;
    }
#ifndef RELEASE
    if( offset != numRecvInd )
        throw std::logic_error("mappedInd was filled incorrectly");
#endif

    // Convert the indices to the original ordering
    inverseMap.Translate( mappedInd );

    // Figure out how many entries each process owns that we need
    mpi::Comm comm = X.Comm();
    const int commSize = mpi::CommSize( comm );
    std::vector<int> recvSizes( commSize, 0 );
    const int blocksize = X.Blocksize();
    for( int s=0; s<numRecvInd; ++s )
    {
        const int i = mappedInd[s];
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
    std::vector<int> recvInd( numRecvInd );
    std::vector<int> offsets = recvOffsets;
    for( int s=0; s<numRecvInd; ++s )
    {
        const int i = mappedInd[s];
        const int q = RowToProcess( i, blocksize, commSize );
        recvInd[offsets[q]++] = i;
    }

    // Coordinate for the coming AllToAll to exchange the indices of X
    std::vector<int> sendSizes( commSize );
    mpi::AllToAll( &recvSizes[0], 1, &sendSizes[0], 1, comm );
    int numSendInd=0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSendInd;
        numSendInd += sendSizes[q];
    }

    // Request the indices
    std::vector<int> sendInd( numSendInd );
    mpi::AllToAll
    ( &recvInd[0], &recvSizes[0], &recvOffsets[0],
      &sendInd[0], &sendSizes[0], &sendOffsets[0], comm );

    // Fulfill the requests
    std::vector<F> sendValues( numSendInd*width_ );
    const int firstLocalRow = X.FirstLocalRow();
    for( int s=0; s<numSendInd; ++s )
        for( int j=0; j<width_; ++j )
            sendValues[s*width_+j] = X.GetLocal( sendInd[s]-firstLocalRow, j );

    // Reply with the values
    std::vector<F> recvValues( numRecvInd*width_ );
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
        const SymmNodeInfo& nodeInfo = info.localNodes[s];
        localNodes[s].ResizeTo( nodeInfo.size, width_ );
        for( int t=0; t<nodeInfo.size; ++t )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width_; ++j )
                localNodes[s].Set( t, j, recvValues[offsets[q]++] );
        }
    }
    distNodes.resize( numDist-1 );
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        DistMatrix<F,VC,STAR>& XNode = distNodes[s-1];
        XNode.SetGrid( *nodeInfo.grid );
        XNode.ResizeTo( nodeInfo.size, width_ );
        const int localHeight = XNode.LocalHeight();
        for( int tLoc=0; tLoc<localHeight; ++tLoc )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width_; ++j )
                XNode.SetLocal( tLoc, j, recvValues[offsets[q]++] );
        }
    }
#ifndef RELEASE
    if( offset != numRecvInd )
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
    const int numSendInd = LocalHeight();
    int offset=0;
    std::vector<int> mappedInd( numSendInd );
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& nodeInfo = info.localNodes[s];
        for( int t=0; t<nodeInfo.size; ++t )
            mappedInd[offset++] = nodeInfo.offset+t;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        const DistMatrix<F,VC,STAR>& XNode = distNodes[s-1];
        for( int t=XNode.ColShift(); t<XNode.Height(); t+=XNode.ColStride() )
            mappedInd[offset++] = nodeInfo.offset+t;
    }

    // Convert the indices to the original ordering
    inverseMap.Translate( mappedInd );

    // Figure out how many indices each process owns that we need to send
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numSendInd; ++s )
    {
        const int i = mappedInd[s];
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
    std::vector<F> sendValues( numSendInd*width );
    std::vector<int> sendInd( numSendInd );
    std::vector<int> offsets = sendOffsets;
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& nodeInfo = info.localNodes[s];
        for( int t=0; t<nodeInfo.size; ++t )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width; ++j )
                sendValues[offsets[q]*width+j] = localNodes[s].Get(t,j);    
            sendInd[offsets[q]++] = i;
        }
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& nodeInfo = info.distNodes[s];
        const DistMatrix<F,VC,STAR>& XNode = distNodes[s-1];
        const int localHeight = XNode.LocalHeight();
        for( int tLoc=0; tLoc<localHeight; ++tLoc )
        {
            const int i = mappedInd[offset++];
            const int q = RowToProcess( i, blocksize, commSize );
            for( int j=0; j<width; ++j )
                sendValues[offsets[q]*width+j] = XNode.GetLocal(tLoc,j);
            sendInd[offsets[q]++] = i;
        }
    }

    // Coordinate for the coming AllToAll to exchange the indices of x
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm );
    int numRecvInd=0;
    std::vector<int> recvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvOffsets[q] = numRecvInd;
        numRecvInd += recvSizes[q];
    }
#ifndef RELEASE
    if( numRecvInd != localHeight )
        throw std::logic_error("numRecvInd was not equal to local height");
#endif

    // Send the indices
    std::vector<int> recvInd( numRecvInd );
    mpi::AllToAll
    ( &sendInd[0], &sendSizes[0], &sendOffsets[0],
      &recvInd[0], &recvSizes[0], &recvOffsets[0], comm );

    // Send the values
    std::vector<F> recvValues( numRecvInd*width );
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
    for( int s=0; s<numRecvInd; ++s )
    {
        const int i = recvInd[s];
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
{
    int localHeight = 0;
    const int numLocal = localNodes.size();
    const int numDist = distNodes.size();
    for( int s=0; s<numLocal; ++s )
        localHeight += localNodes[s].Height();
    for( int s=0; s<numDist; ++s )
        localHeight += distNodes[s].LocalHeight();
    return localHeight;
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_DISTNODALMULTIVEC_IMPL_HPP
