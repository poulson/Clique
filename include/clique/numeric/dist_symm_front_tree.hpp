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
#ifndef CLIQUE_DIST_SYMM_FRONT_TREE_HPP
#define CLIQUE_DIST_SYMM_FRONT_TREE_HPP 1

namespace cliq {

enum SolveMode { NORMAL_1D, NORMAL_2D, FAST_1D_LDL, FAST_2D_LDL };

inline bool 
ModeIs1d( SolveMode mode )
{
    if( mode == NORMAL_1D || mode == FAST_1D_LDL )
        return true;
    else
        return false;
}

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename F>
struct LocalSymmFront
{
    Matrix<F> frontL;
    mutable Matrix<F> work;
};

template<typename F>
struct DistSymmFront
{
    // The 'SolveMode' member variable of the parent 'DistSymmFrontTree' 
    // determines which of the following fronts is active.
    //   {NORMAL_1D,FAST_1D_LDL} -> front1d
    //   {NORMAL_2D,FAST_2D_LDL} -> front2d
    //
    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    // TODO: Think about the fact that almost everything is now mutable...

    mutable DistMatrix<F,VC,STAR> front1dL;
    mutable DistMatrix<F,VC,STAR> work1d;

    mutable DistMatrix<F> front2dL;
    mutable DistMatrix<F> work2d;

    DistMatrix<F,VC,STAR> diag;
};

template<typename F>
struct DistSymmFrontTree
{
    SolveMode mode;
    std::vector<LocalSymmFront<F> > localFronts;
    std::vector<DistSymmFront<F> > distFronts;

    DistSymmFrontTree();

    DistSymmFrontTree
    ( const DistSparseMatrix<F>& A, 
      const std::vector<int>& localMapping,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info,
      Orientation orientation=TRANSPOSE );
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
DistSymmFrontTree<F>::DistSymmFrontTree()
{ }

// TODO: Simplify this using the new MapIndices and InvertMap functions
template<typename F>
DistSymmFrontTree<F>::DistSymmFrontTree
( const DistSparseMatrix<F>& A, 
  const std::vector<int>& localMapping,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info,
  Orientation orientation )
: mode(NORMAL_2D)
{
#ifndef RELEASE
    PushCallStack("DistSymmFrontTree::DistSymmFrontTree");
    if( localMapping.size() != A.NumLocalSources() )
        throw std::logic_error("Local mapping was not the right size");
#endif
    mpi::Comm comm = A.Comm();
    const int blocksize = A.Blocksize();
    const int commSize = mpi::CommSize( comm );

    // Count the number of rows that we want from each process
    std::vector<int> neededRowSizes( commSize, 0 );
    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const int size = info.localNodes[s].size;
        const int offset = info.localNodes[s].offset;
        const std::vector<int>& indices = 
            sepTree.localSepsAndLeaves[s]->indices;
        for( int t=0; t<size; ++t )
        {
            const int i = indices[t];
            const int q = i / blocksize;
            ++neededRowSizes[q];
        }
    }
    const int numDistNodes = info.distNodes.size();
    for( int s=0; s<numDistNodes; ++s )
    {
        // Request entire rows of the distributed matrices, even though we will
        // only need a subset of each row
        
        const DistSymmNodeInfo& node= info.distNodes[s];
        const int size = node.size;
        const int offset = node.offset;
        const int rowShift = node.grid->Col();
        const int rowStride = node.grid->Width();
        const std::vector<int>& indices = sepTree.distSeps[s].indices;
        for( int t=rowShift; t<size; t+=rowStride )
        {
            const int i = indices[t];
            const int q = i / blocksize;
            ++neededRowSizes[q];
        }
    }

    // Fill the set of rows we need from each process
    int numNeededRows=0;
    std::vector<int> neededRowOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        neededRowOffsets[q] = numNeededRows;
        numNeededRows += neededRowSizes[q];
    }
    std::vector<int> neededRows( numNeededRows );
    std::vector<int> offsets = neededRowOffsets;
    for( int s=0; s<numLocalNodes; ++s )
    {
        const int size = info.localNodes[s].size; 
        const int offset = info.localNodes[s].offset;
        const std::vector<int>& indices = 
            sepTree.localSepsAndLeaves[s]->indices;
        for( int t=0; t<size; ++t )
        {
            const int i = indices[t];
            const int q = i / blocksize;
            neededRows[offsets[q]++] = i;
        }
    }
    for( int s=0; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const int size = node.size;
        const int offset = node.offset;
        const int rowShift = node.grid->Col();
        const int rowStride = node.grid->Width();
        const std::vector<int>& indices = sepTree.distSeps[s].indices;
        for( int t=rowShift; t<size; t+=rowStride )
        {
            const int i = indices[t];
            const int q = i / blocksize;
            neededRows[offsets[q]++] = i;
        }
    }

    // Perform an AllToAll to exchange how much each process needs, then 
    // allocate the necessary space for the receives
    std::vector<int> givingRowSizes( commSize );
    mpi::AllToAll( &neededRowSizes[0], 1, &givingRowSizes[0], 1, comm );
    int numGivingRows=0;
    std::vector<int> givingRowOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        givingRowOffsets[q] = numGivingRows;
        numGivingRows += givingRowSizes[q];
    }

    // Perform an AllToAll to exchange what each process needs
    std::vector<int> givingRows( numGivingRows );
    mpi::AllToAll
    ( &neededRows[0], &neededRowSizes[0], &neededRowOffsets[0],
      &givingRows[0], &givingRowSizes[0], &givingRowOffsets[0], comm );

    // Tell each process how many entries we will send per requested row
    std::vector<int> numEntriesGivingPerRow( numGivingRows );
    for( int s=0; s<numGivingRows; ++s )
    {
        const int i = givingRows[s];
        const int firstLocalEntry = A.LocalEntryOffset( i );
        const int lastLocalEntry = A.LocalEntryOffset( i+1 ); 
        numEntriesGivingPerRow[s] = 0;
        for( int t=firstLocalEntry; t<lastLocalEntry; ++t )
            if( A.Col(t) >= i )
                ++numEntriesGivingPerRow[s];
    }
    std::vector<int> numEntriesNeededPerRow( numNeededRows );
    mpi::AllToAll
    ( &numEntriesGivingPerRow[0], &givingRowSizes[0], &givingRowOffsets[0],
      &numEntriesNeededPerRow[0], &neededRowSizes[0], &neededRowOffsets[0],
      comm );

    // Compute the number of entries we give and need from each process
    int numGivingEntries=0;
    std::vector<int> givingEntriesSizes( commSize, 0 );
    std::vector<int> givingEntriesOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = givingRowSizes[q];
        const int offset = givingRowOffsets[q];
        for( int s=0; s<numRows; ++s )
            givingEntriesSizes[q] += numEntriesGivingPerRow[offset+s];

        givingEntriesOffsets[q] = numGivingEntries;
        numGivingEntries += givingEntriesSizes[q];
    }
    numEntriesGivingPerRow.clear();
    int numNeededEntries=0;
    std::vector<int> neededEntriesSizes( commSize, 0 );
    std::vector<int> neededEntriesOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = neededRowSizes[q];
        const int offset = neededRowOffsets[q];
        for( int s=0; s<numRows; ++s )
            neededEntriesSizes[q] += numEntriesNeededPerRow[offset+s];

        neededEntriesOffsets[q] = numNeededEntries;
        numNeededEntries += neededEntriesSizes[q];
    }

    // Pack the information that each process requested
    std::vector<int> givingIndices( numGivingEntries );
    std::vector<F> givingEntries( numGivingEntries );
    offsets = givingEntriesOffsets;
    for( int q=0; q<commSize; ++q )
    {
        int& entryOffset = givingEntriesOffsets[q];
        const int rowOffset = givingRowOffsets[q];
        const int numRows = givingRowSizes[q];
        for( int s=0; s<numRows; ++s )
        {
            const int i = givingRows[rowOffset+s];
            const int firstLocalEntry = A.LocalEntryOffset( i );
            const int lastLocalEntry = A.LocalEntryOffset( i+1 );
            for( int t=firstLocalEntry; t<lastLocalEntry; ++t )
            {
                const int column = A.Col(t);
                if( A.Col(t) >= i )
                {
                    givingIndices[entryOffset] = column;
                    givingEntries[entryOffset] = A.Value(t);
                    ++entryOffset;
                }
            }
        }
    }

    // Perform two AllToAll's to exchange entries and column indices
    std::vector<int> neededIndices( numNeededEntries );
    std::vector<F> neededEntries( numNeededEntries );
    mpi::AllToAll
    ( &givingEntries[0], &givingEntriesSizes[0], &givingEntriesOffsets[0],
      &neededEntries[0], &neededEntriesSizes[0], &neededEntriesOffsets[0], 
      comm );
    givingEntries.clear();
    mpi::AllToAll
    ( &givingIndices[0], &givingEntriesSizes[0], &givingEntriesOffsets[0],
      &neededIndices[0], &neededEntriesSizes[0], &neededEntriesOffsets[0], 
      comm );
    givingIndices.clear();
    givingEntriesSizes.clear();
    givingEntriesOffsets.clear();

    // Perform another round of AllToAll's to get the reordered indices of our
    // recently gathered column indices
    std::set<int> colIndexSet;
    for( int i=0; i<numNeededEntries; ++i )
        colIndexSet.insert( neededIndices[i] );
    const int numNeededMappedIndices = colIndexSet.size();
    std::vector<int> neededMappedIndices( numNeededMappedIndices );
    std::vector<int> neededMappedIndicesOffsets( commSize );
    std::vector<int> neededMappedIndicesSizes( commSize );
    {
        int i=0, qPrev=-1, lastOffset=0;
        std::set<int>::const_iterator it;
        for( it=colIndexSet.begin(); it!=colIndexSet.end(); ++it, ++i )
        {
            const int col = *it;
            neededMappedIndices[i] = col;

            const int q = col / blocksize;
            while( q != qPrev )
            {
                neededMappedIndicesSizes[qPrev] = i - lastOffset;
                neededMappedIndicesOffsets[qPrev+1] = i;
                lastOffset = i;
                ++qPrev;
            }
        }
        neededMappedIndicesSizes[commSize-1] = 
            numNeededMappedIndices-lastOffset;
    }
    std::vector<int> givingMappedIndicesSizes( commSize );
    mpi::AllToAll
    ( &neededMappedIndicesSizes[0], 1, 
      &givingMappedIndicesSizes[0], 1, comm );
    int numGivingMappedIndices=0;
    std::vector<int> givingMappedIndicesOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        givingMappedIndicesOffsets[q] = numGivingMappedIndices;
        numGivingMappedIndices += givingMappedIndicesSizes[q];
    }
    std::vector<int> givingMappedIndices( numGivingMappedIndices );
    mpi::AllToAll
    ( &neededMappedIndices[0], 
      &neededMappedIndicesSizes[0], &neededMappedIndicesOffsets[0],
      &givingMappedIndices[0],
      &givingMappedIndicesSizes[0], &givingMappedIndicesOffsets[0],
      comm );

    // Overwrite givingMappedIndices with their mapped values
    const int firstLocalSource = A.FirstLocalSource();
    for( int i=0; i<numGivingMappedIndices; ++i )
    {
        const int j = givingMappedIndices[i];
        givingMappedIndices[i] = localMapping[j-firstLocalSource];
    }

    // Return the mapped values to everyone
    std::vector<int> mappedIndices( numNeededMappedIndices );
    mpi::AllToAll
    ( &givingMappedIndices[0],
      &givingMappedIndicesSizes[0], &givingMappedIndicesOffsets[0],
      &mappedIndices[0],
      &neededMappedIndicesSizes[0], &neededMappedIndicesOffsets[0],
      comm );
    givingMappedIndices.clear();
    givingMappedIndicesSizes.clear();
    givingMappedIndicesOffsets.clear();

    // Unpack the received information
    const bool conjugate = ( orientation == ADJOINT ); 
    offsets = neededEntriesOffsets;
    std::vector<int> rowOffsets = neededRowOffsets;
    std::vector<int>::const_iterator it;
    localFronts.resize( numLocalNodes );
    for( int s=0; s<numLocalNodes; ++s )
    {
        LocalSymmFront<F>& front = localFronts[s];    
        const LocalSymmNodeInfo& node = info.localNodes[s];

        const int width = node.size;
        const int height = width + node.lowerStruct.size();
        Zeros( height, width, front.frontL );

        const int offset = info.localNodes[s].offset;
        const std::vector<int>& indices = 
            sepTree.localSepsAndLeaves[s]->indices;
        for( int t=0; t<width; ++t )
        {
            const int i = indices[t];
            const int q = i / blocksize;

            const int numEntries = numEntriesNeededPerRow[rowOffsets[q]++];
            for( int k=0; k<numEntries; ++k )
            {
                const int j = neededIndices[offsets[q]];
                it = std::lower_bound
                    ( neededMappedIndices.begin(), 
                      neededMappedIndices.end(), j );
#ifndef RELEASE
                if( it == neededMappedIndices.end() )
                    throw std::logic_error("Did not find needed mapped index");
#endif
                const int indexOffset = it - neededMappedIndices.begin();
                const int mappedIndex = mappedIndices[indexOffset];

                it = std::lower_bound
                    ( node.origLowerStruct.begin(),
                      node.origLowerStruct.end(), mappedIndex );
#ifndef RELEASE
                if( it == node.origLowerStruct.end() )
                    throw std::logic_error("origLowerStruct index not found");
#endif
                const int origOffset = it - node.origLowerStruct.begin();
                const int row = node.origLowerRelIndices[origOffset];

                const F entry = neededEntries[offsets[q]];
                const F value = ( conjugate ? Conj(entry) : entry );
                front.frontL.Set( row, t, value );

                ++offsets[q];
            }
        }
    }
    distFronts.resize( numDistNodes );
    for( int s=0; s<numDistNodes; ++s )
    {
        DistSymmFront<F>& front = distFronts[s];
        const DistSymmNodeInfo& node = info.distNodes[s];
        const Grid& g = *node.grid;

        const int offset = node.offset;
        const int width = node.size;
        const int height = width + node.lowerStruct.size();
        front.front2dL.SetGrid( g );
        Zeros( height, width, front.front2dL );

        const std::vector<int>& indices = sepTree.distSeps[s].indices;
        const int colShift = g.Row();
        const int rowShift = g.Col();
        const int colStride = g.Height();
        const int rowStride = g.Width();
        for( int t=rowShift; t<width; t+=rowStride )
        {
            const int i = indices[t];
            const int q = i / blocksize;

            const int numEntries = numEntriesNeededPerRow[rowOffsets[q]++];
            for( int k=0; k<numEntries; ++k )
            {
                const int j = neededIndices[offsets[q]];
                it = std::lower_bound
                    ( neededMappedIndices.begin(), 
                      neededMappedIndices.end(), j );
#ifndef RELEASE
                if( it == neededMappedIndices.end() )
                    throw std::logic_error("Did not find needed mapped index");
#endif
                const int indexOffset = it - neededMappedIndices.begin();
                const int mappedIndex = mappedIndices[indexOffset];

                it = std::lower_bound
                    ( node.origLowerStruct.begin(), 
                      node.origLowerStruct.end(), mappedIndex );
#ifndef RELEASE
                if( it == node.origLowerStruct.end() )
                    throw std::logic_error("origLowerStruct index not found");
#endif
                const int origOffset = it - node.origLowerStruct.begin();
                const int row = node.origLowerRelIndices[origOffset];

                if( row % colStride == colShift )
                {
                    const int localRow = (row-colShift) / colStride;
                    const int localCol = (t-rowShift) / rowStride;
                    const F entry = neededEntries[offsets[q]];
                    const F value = ( conjugate ? Conj(entry) : entry );
                    front.frontL.SetLocal( localRow, localCol, value );
                }

                ++offsets[q];
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // CLIQUE_DIST_SYMM_FRONT_TREE_HPP
