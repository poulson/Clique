/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CLIQUE_DIST_NODAL_VECTOR_HPP
#define CLIQUE_DIST_NODAL_VECTOR_HPP 1

namespace cliq {

// For handling a vector distributed in a [VC,* ] manner over each node
// of the elimination tree
template<typename F>
class DistNodalVector
{
public:
    std::vector<F> values;

    DistNodalVector
    ( const std::vector<int>& localInverseMap, const DistSymmInfo& info, 
      const DistVector<F>& x );
private:
    // This is made private to prevent default constructors from being called
    DistNodalVector();
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
DistNodalVector<F>::DistNodalVector
( const std::vector<int>& localInverseMap, const DistSymmInfo& info,
  const DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("DistNodalVector::DistNodalVector");
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
        numRecvIndices += node.size;
    }
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
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

    // Convert the indices to the original ordering
    MapIndices( localInverseMap, mappedIndices, x.Height(), comm );

    // Figure out how many entries each process owns that we need
    std::vector<int> recvIndexSizes( commSize );
    for( int s=0; s<numRecvIndices; ++s )
    {
        const int i = mappedIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++recvIndexSizes[q];
    }
    std::vector<int> recvIndexOffsets( commSize );
    {
        offset=0;
        for( int q=0; q<commSize; ++q )
        {
            recvIndexOffsets[q] = offset;
            offset += recvIndexSizes[q];
        }
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
    values.resize( numRecvIndices );
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        for( int t=0; t<node.size; ++t )
        {
            const int i = mappedIndices[offset];
            const int q = RowToProcess( i, blocksize, commSize );
            values[offset++] = recvValues[offsets[q]++];
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
            values[offset++] = recvValues[offsets[q]++];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // CLIQUE_DIST_NODAL_VECTOR_HPP
