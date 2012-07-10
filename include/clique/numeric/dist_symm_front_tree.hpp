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

template<typename F>
DistSymmFrontTree<F>::DistSymmFrontTree
( const DistSparseMatrix<F>& A, 
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info,
  Orientation orientation )
: mode(NORMAL_2D)
{
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    const int numLocalNodes = info.localNodes.size();
    const int numDistNodes = info.distNodes.size();
    
    localFronts.resize( numLocalNodes );
    distFronts.resize( numDistNodes );

    //
    // Generate the description of the data that we want from each process
    //
    const int blocksize = A.Blocksize();

    // Count the number of rows that we want from each process
    std::vector<int> numNeededRows( commSize, 0 );
    for( int s=0; s<numLocalNodes; ++s )
    {
        const int size = info.localNodes[s].size;
        const int offset = info.localNodes[s].offset;
        const std::vector<int>& indices = sepTree.localSepsAndLeaves[s].indices;
        for( int t=0; t<size; ++t )
        {
            const int i = indices[t];
            const int proc = i / blocksize;
            ++numNeededRows[proc];
        }
    }
    for( int s=0; s<numDistNodes; ++s )
    {
        // Request entire rows of the distributed matrices, even though we will
        // only need a subset of them
        const int size = info.distNodes[s].size; 
        const int offset = info.distNodes[s].offset;
        const std::vector<int>& indices = sepTree.distSeps[s].indices;
        for( int t=0; t<size; ++t )
        {
            const int i = indices[t];
            const int proc = i / blocksize;
            ++numNeededRows[proc];
        }
    }
    // Fill the set of rows we need from each process
    int totalNumNeededRows=0;
    std::vector<int> neededRowOffsets( commSize );
    for( int i=0; i<commSize; ++i )
    {
        neededRowOffsets[i] = totalNumNeededRows;
        totalNumNeededRows += numNeededRows[i];
    }
    std::vector<int> neededRows( totalNumNeededRows );
    std::vector<int> offsets = neededRowOffsets;
    for( int s=0; s<numLocalNodes; ++s )
    {
        const int size = info.localNodes[s].size; 
        const int offset = info.localNodes[s].offset;
        const std::vector<int>& indices = sepTree.localSepsAndLeaves[s].indices;
        for( int t=0; t<size; ++t )
        {
            const int i = indices[t];
            const int proc = i / blocksize;
            neededRows[offsets[proc]] = i;
            ++offsets[proc];
        }
    }
    for( int s=0; s<numDistNodes; ++s )
    {
        const int size = info.distNodes[s].size;
        const int offset = info.distNodes[s].offset;
        const std::vector<int>& indices = sepTree.distSeps[s].indices;
        for( int t=0; t<size; ++t )
        {
            const int i = indices[t];
            const int proc = i / blocksize;
            neededRows[offsets[proc]] = i;
            ++offsets[proc];
        }
    }

    // Perform an AllToAll to exchange how much each process needs, then 
    // allocate the necessary space for the receives
    std::vector<int> numProvidedRows( commSize );
    mpi::AllToAll( &numNeededRows[0], 1, &numProvidedRows[0], 1, comm );
    int totalNumProvidedRows=0;
    std::vector<int> providedRowOffsets( commSize );
    for( int i=0; i<commSize; ++i )
    {
        providedRowOffsets[i] = totalNumProvidedRows;
        totalNumProvidedRows += numProvidedRows[i];
    }
    std::vector<int> providedRows( totalNumProvidedRows );

    // Perform an AllToAll to exchange what each process needs
    mpi::AllToAll
    ( &neededRows[0], &numNeededRows[0], &neededRowOffsets[0],
      &providedRows[0], &numProvidedRows[0], &providedRowOffsets[0], comm );

    // Tell each process how many entries we will send per requested row
    std::vector<int> numEntriesProvidedPerRow( totalNumProvidedRows );
    for( int s=0; s<totalNumProvidedRows; ++s )
    {
        const int i = providedRows[s];
        const int firstLocalEntry = A.LocalEntryOffset( i );
        const int lastLocalEntry = A.LocalEntryOffset( i+1 ); 
        numEntriesProvidedPerRow[s] = 0;
        for( int t=firstLocalEntry; t<lastLocalEntry; ++t )
            if( A.Col(t) >= i )
                ++numEntriesProvidedPerRow[s];
    }
    std::vector<int> numEntriesNeededPerRow( totalNumNeededRows );
    mpi::AllToAll
    ( &numEntriesProvidedPerRow[0], &numProvidedRows[0], &providedRowOffsets[0],
      &numEntriesNeededPerRow[0], &numNeededRows[0], &neededRowOffsets[0],
      comm );

    // Pack the information that each process requested
    // TODO

    // Perform an AllToAll to exchange what each process needed
    // TODO

    // Unpack the received information
    const bool conjugate = ( orientation == ADJOINT );
    // TODO
    for( int s=0; s<numLocalNodes; ++s )
    {
        LocalSymmFront<F>& front = localFronts[s];    
        const LocalSymmNodeInfo& frontInfo = info.localNodes[s];

        const int width = frontInfo.size;
        const int height = width + frontInfo.lowerStruct.size();
        Zeros( height, width, front.frontL );
        // TODO: Add onto front
    }
    for( int s=0; s<numDistNodes; ++s )
    {
        DistSymmFront<F>& front = distFronts[s];
        const DistSymmNodeInfo& frontInfo = info.distNodes[s];

        const int width = frontInfo.size;
        const int height = width + frontInfo.lowerStruct.size();
        front.front2dL.SetGrid( *frontInfo.grid );
        Zeros( height, width, front.front2dL );
        // TODO: Add onto front
    }
}

} // namespace cliq

#endif // CLIQUE_DIST_SYMM_FRONT_TREE_HPP
