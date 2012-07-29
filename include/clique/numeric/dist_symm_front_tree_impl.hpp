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
DistSymmFrontTree<F>::DistSymmFrontTree()
{ }

template<typename F>
DistSymmFrontTree<F>::DistSymmFrontTree
( Orientation orientation,
  const DistSparseMatrix<F>& A, 
  const std::vector<int>& localMap,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info )
: mode(NORMAL_2D)
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::DistSymmFrontTree");
    if( A.LocalHeight() != (int)localMap.size() )
        throw std::logic_error("Local mapping was not the right size");
#endif
    mpi::Comm comm = A.Comm();
    const DistGraph& graph = A.Graph();
    const int blocksize = A.Blocksize();
    const int commSize = mpi::CommSize( comm );
    const int numSources = graph.NumSources();
    const int numLocal = sepTree.localSepsAndLeaves.size();
    const int numDist = sepTree.distSeps.size();

    // Get the reordered indices of the targets of our portion of the 
    // distributed sparse matrix
    std::set<int> targetSet( graph.targets_.begin(), graph.targets_.end() );
    std::vector<int> targets( targetSet.size() );
    std::copy( targetSet.begin(), targetSet.end(), targets.begin() );
    std::vector<int> mappedTargets = targets;
    MapIndices( localMap, mappedTargets, numSources, comm );

    // Set up the indices for the rows we need from each process
    std::vector<int> recvRowSizes( commSize, 0 );
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const int numIndices = sepOrLeaf.indices.size();
        for( int t=0; t<numIndices; ++t )
        {
            const int i = sepOrLeaf.indices[t];
#ifndef RELEASE
            if( i < 0 || i >= numSources )
                throw std::logic_error("separator index was out of bounds");
#endif
            const int q = RowToProcess( i, blocksize, commSize );
            ++recvRowSizes[q];
        }
    }
    for( int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const Grid& grid = *node.grid;
        const int rowShift = grid.Col();
        const int rowStride = grid.Width();
        const int numIndices = sep.indices.size();
        for( int t=rowShift; t<numIndices; t+=rowStride )
        {
            const int i = sep.indices[t];
#ifndef RELEASE
            if( i < 0 || i >= numSources )
                throw std::logic_error("separator index was out of bounds");
#endif
            const int q = RowToProcess( i, blocksize, commSize );
            ++recvRowSizes[q];
        }
    }
    int numRecvRows=0;
    std::vector<int> recvRowOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvRowOffsets[q] = numRecvRows;
        numRecvRows += recvRowSizes[q];
    }
    std::vector<int> recvRows( numRecvRows );
    std::vector<int> offsets = recvRowOffsets;
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const int numIndices = sepOrLeaf.indices.size();
        for( int t=0; t<numIndices; ++t )
        {
            const int i = sepOrLeaf.indices[t];
#ifndef RELEASE
            if( i < 0 || i >= numSources )
                throw std::logic_error("separator index was out of bounds");
#endif
            const int q = RowToProcess( i, blocksize, commSize );
#ifndef RELEASE            
            if( offsets[q] >= numRecvRows )
                throw std::logic_error("offset got too large");
#endif
            recvRows[offsets[q]++] = i;
        }
    }
    for( int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const Grid& grid = *node.grid;
        const int rowShift = grid.Col();
        const int rowStride = grid.Width();
        const int numIndices = sep.indices.size();
        for( int t=rowShift; t<numIndices; t+=rowStride )
        {
            const int i = sep.indices[t];
#ifndef RELEASE
            if( i < 0 || i >= numSources )
                throw std::logic_error("separator index was out of bounds");
#endif
            const int q = RowToProcess( i, blocksize, commSize );
#ifndef RELEASE            
            if( offsets[q] >= numRecvRows )
                throw std::logic_error("offset got too large");
#endif
            recvRows[offsets[q]++] = i;
        }
    }

    // Retreive the list of rows that we must send to each process
    std::vector<int> sendRowSizes( commSize );
    mpi::AllToAll( &recvRowSizes[0], 1, &sendRowSizes[0], 1, comm );
    int numSendRows=0;
    std::vector<int> sendRowOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendRowOffsets[q] = numSendRows;
        numSendRows += sendRowSizes[q];
    }
    std::vector<int> sendRows( numSendRows );
    mpi::AllToAll
    ( &recvRows[0], &recvRowSizes[0], &recvRowOffsets[0],
      &sendRows[0], &sendRowSizes[0], &sendRowOffsets[0], comm );

    // Pack the number of nonzeros per row (and the nonzeros themselves)
    int numSendEntries=0;
    const int firstLocalRow = A.FirstLocalRow();
    std::vector<int> sendRowLengths( numSendRows );
    std::vector<int> sendEntriesSizes( commSize, 0 );
    std::vector<int> sendEntriesOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int size = sendRowSizes[q];
        const int offset = sendRowOffsets[q];
        sendEntriesOffsets[q] = numSendEntries;
        for( int s=0; s<size; ++s )
        {
            const int i = sendRows[s+offset];
            const int iLocal = i - firstLocalRow;
            const int numConnections = A.NumConnections( iLocal );
            numSendEntries += numConnections;
            sendEntriesSizes[q] += numConnections;
            sendRowLengths[s+offset] = numConnections;
        }
    }
    const bool conjugate = ( orientation == ADJOINT ? true : false );
    std::vector<F> sendEntries( numSendEntries );
    std::vector<int> sendTargets( numSendEntries );
    std::vector<int>::const_iterator vecIt;
    for( int q=0; q<commSize; ++q )
    {
        int index = sendEntriesOffsets[q];
        const int size = sendRowSizes[q];
        const int offset = sendRowOffsets[q];
        for( int s=0; s<size; ++s )
        {
            const int i = sendRows[s+offset];
            const int iLocal = i - firstLocalRow;
            const int numConnections = sendRowLengths[s+offset];
            const int localEntryOffset = A.LocalEntryOffset( iLocal );
            for( int t=0; t<numConnections; ++t )
            {
                const F value = A.Value( localEntryOffset+t );
                const int col = A.Col( localEntryOffset+t );
                vecIt = std::lower_bound
                    ( targets.begin(), targets.end(), col );
#ifndef RELEASE
                if( vecIt == targets.end() )
                    throw std::logic_error("Could not find target");
#endif
                const int targetOffset = vecIt - targets.begin();
                const int mappedTarget = mappedTargets[targetOffset];
#ifndef RELEASE
                if( index >= numSendEntries )
                    throw std::logic_error("send entry index got too big");
#endif
                sendEntries[index] = ( conjugate ? elem::Conj(value) : value );
                sendTargets[index] = mappedTarget;
                ++index;
            }
        }
#ifndef RELEASE
        if( index != sendEntriesOffsets[q]+sendEntriesSizes[q] )
            throw std::logic_error("index was not the correct value");
#endif
    }

    // Send back the number of nonzeros per row and the nonzeros themselves
    std::vector<int> recvRowLengths( numRecvRows );
    mpi::AllToAll
    ( &sendRowLengths[0], &sendRowSizes[0], &sendRowOffsets[0],
      &recvRowLengths[0], &recvRowSizes[0], &recvRowOffsets[0], comm );
    int numRecvEntries=0;
    std::vector<int> recvEntriesSizes( commSize, 0 );
    std::vector<int> recvEntriesOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int size = recvRowSizes[q];
        const int offset = recvRowOffsets[q];
        for( int s=0; s<size; ++s )
            recvEntriesSizes[q] += recvRowLengths[offset+s];

        recvEntriesOffsets[q] = numRecvEntries; 
        numRecvEntries += recvEntriesSizes[q];
    }
    std::vector<F> recvEntries( numRecvEntries );
    std::vector<int> recvTargets( numRecvEntries );
    mpi::AllToAll
    ( &sendEntries[0], &sendEntriesSizes[0], &sendEntriesOffsets[0],
      &recvEntries[0], &recvEntriesSizes[0], &recvEntriesOffsets[0], comm );
    mpi::AllToAll
    ( &sendTargets[0], &sendEntriesSizes[0], &sendEntriesOffsets[0],
      &recvTargets[0], &recvEntriesSizes[0], &recvEntriesOffsets[0], comm );

    // Unpack the received entries
    offsets = recvRowOffsets;
    std::vector<int> entryOffsets = recvEntriesOffsets;
    localFronts.resize( numLocal );
    for( int s=0; s<numLocal; ++s )
    {
        LocalSymmFront<F>& front = localFronts[s];
        const LocalSepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const LocalSymmNodeInfo& node = info.localNodes[s];
        const std::vector<int>& origLowerStruct = node.origLowerStruct;

        const int size = node.size;
        const int offset = node.offset;
        const int lowerSize = node.lowerStruct.size();
        Zeros( size+lowerSize, size, front.frontL );

#ifndef RELEASE
        if( size != (int)sepOrLeaf.indices.size() )
            throw std::logic_error("Mismatch between separator and node size");
#endif

        for( int t=0; t<size; ++t )
        {
            const int i = sepOrLeaf.indices[t];
            const int q = RowToProcess( i, blocksize, commSize );

            int& entryOffset = entryOffsets[q];
            const int numEntries = recvRowLengths[offsets[q]++];

            for( int k=0; k<numEntries; ++k )
            {
                const F value = recvEntries[entryOffset];
                const int target = recvTargets[entryOffset];
                ++entryOffset;

                if( target < offset+t )
                    continue;
                else if( target < offset+size )
                {
                    front.frontL.Set( target-offset, t, value );
                }
                else
                {
                    vecIt = std::lower_bound
                        ( origLowerStruct.begin(), 
                          origLowerStruct.end(), target );
#ifndef RELEASE
                    if( vecIt == origLowerStruct.end() )
                        throw std::logic_error("No match in origLowerStruct");
#endif
                    const int origOffset = vecIt - origLowerStruct.begin();
#ifndef RELEASE
                    if( origOffset >= (int)node.origLowerRelIndices.size() )
                        throw std::logic_error("origLowerRelIndices too small");
#endif
                    const int row = node.origLowerRelIndices[origOffset];
                    front.frontL.Set( row, t, value );
                }
            }
        }
    }

    distFronts.resize( numDist+1 );
    for( int s=0; s<numDist; ++s )
    {
        DistSymmFront<F>& front = distFronts[s+1];
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const std::vector<int>& origLowerStruct = node.origLowerStruct;

        const Grid& grid = *node.grid;
        const int colShift = grid.Row();
        const int rowShift = grid.Col();
        const int colStride = grid.Height();
        const int rowStride = grid.Width();

        const int size = node.size;
        const int offset = node.offset;
        const int lowerSize = node.lowerStruct.size();
        front.front2dL.SetGrid( grid );
        Zeros( size+lowerSize, size, front.front2dL );

#ifndef RELEASE
        if( size != (int)sep.indices.size() )
            throw std::logic_error("Mismatch in separator and node sizes");
#endif

        for( int t=rowShift; t<size; t+=rowStride )
        {
            const int i = sep.indices[t];
            const int q = RowToProcess( i, blocksize, commSize );
            const int localCol = (t-rowShift) / rowStride;

            int& entryOffset = entryOffsets[q];
            const int numEntries = recvRowLengths[offsets[q]++];

            for( int k=0; k<numEntries; ++k )
            {
                const F value = recvEntries[entryOffset];
                const int target = recvTargets[entryOffset];
                ++entryOffset;

                if( target < offset+t )
                    continue;
                else if( target < offset+size )
                {
                    if( (target-offset) % colStride == colShift )
                    {
                        const int row = target-offset;
                        const int localRow = (row-colShift) / colStride;
                        front.front2dL.SetLocal( localRow, localCol, value );
                    }
                }
                else 
                {
                    vecIt = std::lower_bound
                        ( origLowerStruct.begin(),
                          origLowerStruct.end(), target );
#ifndef RELEASE
                    if( vecIt == origLowerStruct.end() )
                        throw std::logic_error("No match in origLowerStruct");
#endif
                    const int origOffset = vecIt - origLowerStruct.begin();
#ifndef RELEASE
                    if( origOffset >= (int)node.origLowerRelIndices.size() )
                        throw std::logic_error("origLowerRelIndices too small");
#endif
                    const int row = node.origLowerRelIndices[origOffset];
                    if( row % colStride == colShift )
                    {
                        const int localRow = (row-colShift) / colStride;
                        front.front2dL.SetLocal( localRow, localCol, value );
                    }
                }
            }
        }
    }
#ifndef RELEASE
    for( int q=0; q<commSize; ++q )
        if( entryOffsets[q] != recvEntriesOffsets[q]+recvEntriesSizes[q] )
            throw std::logic_error("entryOffsets were incorrect");
#endif
    
    // Copy information from the local root to the dist leaf
    {
        const DistSymmNodeInfo& node = info.distNodes[0];
        Matrix<F>& topLocal = localFronts.back().frontL;
        DistMatrix<F>& bottomDist = distFronts[0].front2dL;
        bottomDist.LockedView
        ( topLocal.Height(), topLocal.Width(), 0, 0,
          topLocal.LockedBuffer(), topLocal.LDim(), *node.grid );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
