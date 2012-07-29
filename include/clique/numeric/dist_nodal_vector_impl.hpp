/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
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

namespace cliq {

template<typename F>
inline
DistNodalVector<F>::DistNodalVector()
{ }

template<typename F>
inline
DistNodalVector<F>::DistNodalVector
( const std::vector<int>& localInverseMap, const DistSymmInfo& info,
  const DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("DistNodalVector::DistNodalVector");
#endif
    Pull( localInverseMap, info, x );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistNodalVector<F>::Pull
( const std::vector<int>& localInverseMap, const DistSymmInfo& info,
  const DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("DistNodalVector::Pull");
#endif
    mpi::Comm comm = x.Comm();
    const int commSize = mpi::CommSize( comm );
    const int blocksize = x.Blocksize();
    const int firstLocalRow = x.FirstLocalRow();
    const int numDist = info.distNodes.size();
    const int numLocal = info.localNodes.size();

    // Traverse our part of the elimination tree to see how many indices we need
    int numRecvIndices=0;
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
#ifndef RELEASE
        if( numRecvIndices != node.myOffset )
            throw std::logic_error("numRecvIndices did not match local offset");
#endif
        numRecvIndices += node.size;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
#ifndef RELEASE
        if( numRecvIndices != node.localOffset1d )
            throw std::logic_error("numRecvIndices did not match dist offset");
#endif
        numRecvIndices += node.localSize1d;
    }
    
    // Fill the set of indices that we need to map to the original ordering
    int offset=0;
    std::vector<int> mappedIndices( numRecvIndices );
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
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
    MapIndices( localInverseMap, mappedIndices, x.Height(), comm );

    // Figure out how many entries each process owns that we need
    std::vector<int> recvIndexSizes( commSize, 0 );
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++recvIndexSizes[q];
    }
    std::vector<int> recvIndexOffsets( commSize );
    offset=0;
    for( int q=0; q<commSize; ++q )
    {
        recvIndexOffsets[q] = offset;
        offset += recvIndexSizes[q];
    }
    std::vector<int> recvIndices( numRecvIndices );
    std::vector<int> offsets = recvIndexOffsets;
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        recvIndices[offsets[q]++] = i;
    }

    // Coordinate for the coming AllToAll to exchange the indices of x
    std::vector<int> sendIndexSizes( commSize );
    mpi::AllToAll( &recvIndexSizes[0], 1, &sendIndexSizes[0], 1, comm );
    int numSendIndices=0;
    std::vector<int> sendIndexOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendIndexOffsets[q] = numSendIndices;
        numSendIndices += sendIndexSizes[q];
    }

    // Request the indices
    std::vector<int> sendIndices( numSendIndices );
    mpi::AllToAll
    ( &recvIndices[0], &recvIndexSizes[0], &recvIndexOffsets[0],
      &sendIndices[0], &sendIndexSizes[0], &sendIndexOffsets[0], comm );

    // Fulfill the requests
    std::vector<F> sendValues( numSendIndices );
    for( int s=0; s<numSendIndices; ++s )
        sendValues[s] = x.GetLocal( sendIndices[s]-firstLocalRow );

    // Reply with the values
    std::vector<F> recvValues( numRecvIndices );
    mpi::AllToAll
    ( &sendValues[0], &sendIndexSizes[0], &sendIndexOffsets[0],
      &recvValues[0], &recvIndexSizes[0], &recvIndexOffsets[0], comm );
    sendValues.clear();
    sendIndexSizes.clear();
    sendIndexOffsets.clear();

    // Unpack the values
    offset=0;
    offsets = recvIndexOffsets;
    localVec.ResizeTo( numRecvIndices, 1 );
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
        {
            const int i = mappedIndices[offset];
            const int q = RowToProcess( i, blocksize, commSize );
#ifndef RELEASE
            if( offsets[q] >= numRecvIndices )
                throw std::logic_error("One of the offsets was too large");
#endif
            localVec.Set( offset++, 0, recvValues[offsets[q]++] );
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
        for( int t=shift; t<node.size; t+=gridSize )
        {
            const int i = mappedIndices[offset];
            const int q = RowToProcess( i, blocksize, commSize );
#ifndef RELEASE
            if( offsets[q] >= numRecvIndices )
                throw std::logic_error("One of the offsets was too large");
#endif
            localVec.Set( offset++, 0, recvValues[offsets[q]++] );
        }
    }
#ifndef RELEASE
    if( offset != numRecvIndices )
        throw std::logic_error("Unpacked wrong number of indices");
#endif

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistNodalVector<F>::Push
( const std::vector<int>& localInverseMap, const DistSymmInfo& info,
        DistVector<F>& x ) const
{
#ifndef RELEASE
    PushCallStack("DistNodalVector::Push");
#endif
    const DistSymmNodeInfo& rootNode = info.distNodes.back();
    mpi::Comm comm = rootNode.comm;
    const int height = rootNode.size + rootNode.offset;
    x.SetComm( comm );
    x.ResizeTo( height );

    const int commSize = mpi::CommSize( comm );
    const int blocksize = x.Blocksize();
    const int localHeight = x.LocalHeight();
    const int firstLocalRow = x.FirstLocalRow();
    const int numDist = info.distNodes.size();
    const int numLocal = info.localNodes.size();

    // Fill the set of indices that we need to map to the original ordering
    const int numSendIndices = localVec.Height();
    int offset=0;
    std::vector<int> mappedIndices( numSendIndices );
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
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
    MapIndices( localInverseMap, mappedIndices, x.Height(), comm );

    // Figure out how many entries each process owns that we need to send
    std::vector<int> sendIndexSizes( commSize, 0 );
    for( int s=0; s<numSendIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++sendIndexSizes[q];
    }
    std::vector<int> sendIndexOffsets( commSize );
    offset=0;
    for( int q=0; q<commSize; ++q )
    {
        sendIndexOffsets[q] = offset;
        offset += sendIndexSizes[q];
    }

    // Pack the send indices and values
    offset=0;
    std::vector<F> sendValues( numSendIndices );
    std::vector<int> sendIndices( numSendIndices );
    std::vector<int> offsets = sendIndexOffsets;
    for( int s=0; s<numSendIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        sendValues[offsets[q]] = localVec.Get(offset++,0);
        sendIndices[offsets[q]] = i;
        ++offsets[q];
    }

    // Coordinate for the coming AllToAll to exchange the indices of x
    std::vector<int> recvIndexSizes( commSize );
    mpi::AllToAll( &sendIndexSizes[0], 1, &recvIndexSizes[0], 1, comm );
    int numRecvIndices=0;
    std::vector<int> recvIndexOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvIndexOffsets[q] = numRecvIndices;
        numRecvIndices += recvIndexSizes[q];
    }
#ifndef RELEASE
    if( numRecvIndices != localHeight )
        throw std::logic_error("numRecvIndices was not equal to local height");
#endif

    // Send the indices
    std::vector<int> recvIndices( numRecvIndices );
    mpi::AllToAll
    ( &sendIndices[0], &sendIndexSizes[0], &sendIndexOffsets[0],
      &recvIndices[0], &recvIndexSizes[0], &recvIndexOffsets[0], comm );

    // Send the values
    std::vector<F> recvValues( numRecvIndices );
    mpi::AllToAll
    ( &sendValues[0], &sendIndexSizes[0], &sendIndexOffsets[0],
      &recvValues[0], &recvIndexSizes[0], &recvIndexOffsets[0], comm );
    sendValues.clear();
    sendIndexSizes.clear();
    sendIndexOffsets.clear();

    // Unpack the values
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = recvIndices[s];
        const int iLocal = i - firstLocalRow;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= localHeight )
            throw std::logic_error("iLocal was out of bounds");
#endif
        x.SetLocal( iLocal, recvValues[s] );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
