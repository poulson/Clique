/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F>
inline
DistSymmFrontTree<F>::DistSymmFrontTree()
{ }

template<typename F>
inline void
DistSymmFrontTree<F>::Initialize
( Orientation orientation,
  const DistSparseMatrix<F>& A, 
  const DistMap& reordering,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info )
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::Initialize");
    if( orientation == NORMAL )
        throw std::logic_error("Matrix must be symmetric or Hermitian");
    if( A.LocalHeight() != reordering.NumLocalSources() )
        throw std::logic_error("Local mapping was not the right size");
#endif
    isHermitian = ( orientation != TRANSPOSE );
    frontType = SYMM_2D;
    
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
    reordering.Translate( mappedTargets );

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
    std::vector<F> sendEntries( numSendEntries );
    std::vector<int> sendTargets( numSendEntries );
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
                const int targetOffset = Find( targets, col );
                const int mappedTarget = mappedTargets[targetOffset];
#ifndef RELEASE
                if( index >= numSendEntries )
                    throw std::logic_error("send entry index got too big");
#endif
                sendEntries[index] = (isHermitian ? elem::Conj(value) : value);
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
                    const int origOffset = Find( origLowerStruct, target );
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
                    const int origOffset = Find( origLowerStruct, target );
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

template<typename F>
inline
DistSymmFrontTree<F>::DistSymmFrontTree
( Orientation orientation,
  const DistSparseMatrix<F>& A, 
  const DistMap& reordering,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info )
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::DistSymmFrontTree");
#endif
    Initialize( orientation, A, reordering, sepTree, info );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistSymmFrontTree<F>::MemoryInfo
( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries, 
  double& numGlobalEntries ) const
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::MemInfo");
#endif
    numLocalEntries = numGlobalEntries = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const LocalSymmFront<F>& front = localFronts[s];
        numLocalEntries += front.frontL.MemorySize();
        numLocalEntries += front.work.MemorySize();
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<F>& front = distFronts[s];
        if( frontsAre1d )
        {
            numLocalEntries += front.front1dL.AllocatedMemory();
            numLocalEntries += front.work1d.AllocatedMemory();
        }
        else
        {
            numLocalEntries += front.front2dL.AllocatedMemory();
            numLocalEntries += front.work2d.AllocatedMemory();
        }
        numLocalEntries += front.diag.AllocatedMemory();
    }

    mpi::AllReduce( &numLocalEntries, &minLocalEntries, 1, mpi::MIN, comm );
    mpi::AllReduce( &numLocalEntries, &maxLocalEntries, 1, mpi::MAX, comm );
    mpi::AllReduce( &numLocalEntries, &numGlobalEntries, 1, mpi::SUM, comm );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistSymmFrontTree<F>::TopLeftMemoryInfo
( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries, 
  double& numGlobalEntries ) const
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::TopLeftMemInfo");
#endif
    numLocalEntries = numGlobalEntries = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const LocalSymmFront<F>& front = localFronts[s];
        Matrix<F> FTL,
                  FBL;
        elem::LockedPartitionDown
        ( front.frontL, FTL,
                        FBL, front.frontL.Width() );
        numLocalEntries += FTL.Height()*FTL.Width();
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<F>& front = distFronts[s];
        if( frontsAre1d )
        {
            DistMatrix<F,VC,STAR> FTL(grid),
                                  FBL(grid);
            elem::LockedPartitionDown
            ( front.front1dL, FTL,
                              FBL, front.front1dL.Width() );
            numLocalEntries += FTL.LocalHeight()*FTL.LocalWidth();
        }
        else
        {
            DistMatrix<F> FTL(grid),
                          FBL(grid);
            elem::LockedPartitionDown
            ( front.front2dL, FTL,
                              FBL, front.front2dL.Width() );
            numLocalEntries += FTL.LocalHeight()*FTL.LocalWidth();
        }
        numLocalEntries += front.diag.AllocatedMemory();
    }

    mpi::AllReduce( &numLocalEntries, &minLocalEntries, 1, mpi::MIN, comm );
    mpi::AllReduce( &numLocalEntries, &maxLocalEntries, 1, mpi::MAX, comm );
    mpi::AllReduce( &numLocalEntries, &numGlobalEntries, 1, mpi::SUM, comm );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistSymmFrontTree<F>::BottomLeftMemoryInfo
( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries, 
  double& numGlobalEntries ) const
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::BottomLeftMemInfo");
#endif
    numLocalEntries = numGlobalEntries = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const LocalSymmFront<F>& front = localFronts[s];
        Matrix<F> FTL,
                  FBL;
        elem::LockedPartitionDown
        ( front.frontL, FTL,
                        FBL, front.frontL.Width() );
        numLocalEntries += FBL.Height()*FBL.Width();
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<F>& front = distFronts[s];
        if( frontsAre1d )
        {
            DistMatrix<F,VC,STAR> FTL(grid),
                                  FBL(grid);
            elem::LockedPartitionDown
            ( front.front1dL, FTL,
                              FBL, front.front1dL.Width() );
            numLocalEntries += FBL.LocalHeight()*FBL.LocalWidth();
        }
        else
        {
            DistMatrix<F> FTL(grid),
                          FBL(grid);
            elem::LockedPartitionDown
            ( front.front2dL, FTL,
                              FBL, front.front2dL.Width() );
            numLocalEntries += FBL.LocalHeight()*FBL.LocalWidth();
        }
        numLocalEntries += front.diag.AllocatedMemory();
    }

    mpi::AllReduce( &numLocalEntries, &minLocalEntries, 1, mpi::MIN, comm );
    mpi::AllReduce( &numLocalEntries, &maxLocalEntries, 1, mpi::MAX, comm );
    mpi::AllReduce( &numLocalEntries, &numGlobalEntries, 1, mpi::SUM, comm );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
