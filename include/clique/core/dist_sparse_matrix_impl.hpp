/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T>
inline 
DistSparseMatrix<T>::DistSparseMatrix()
{ }

template<typename T>
inline 
DistSparseMatrix<T>::DistSparseMatrix( mpi::Comm comm )
: graph_(comm)
{ }

template<typename T>
inline
DistSparseMatrix<T>::DistSparseMatrix( int height, mpi::Comm comm )
: graph_(height,comm)
{ }

template<typename T>
inline 
DistSparseMatrix<T>::DistSparseMatrix( int height, int width, mpi::Comm comm )
: graph_(height,width,comm)
{ }

template<typename T>
inline 
DistSparseMatrix<T>::~DistSparseMatrix()
{ }

template<typename T>
inline int 
DistSparseMatrix<T>::Height() const
{ return graph_.NumSources(); }

template<typename T>
inline int 
DistSparseMatrix<T>::Width() const
{ return graph_.NumTargets(); }

template<typename T>
inline const DistGraph& 
DistSparseMatrix<T>::Graph() const
{ return graph_; }

template<typename T>
inline void
DistSparseMatrix<T>::SetComm( mpi::Comm comm )
{ 
    graph_.SetComm( comm ); 
    values_.clear();
}

template<typename T>
inline mpi::Comm 
DistSparseMatrix<T>::Comm() const
{ return graph_.Comm(); }

template<typename T>
inline int
DistSparseMatrix<T>::Blocksize() const
{ return graph_.Blocksize(); }

template<typename T>
inline int
DistSparseMatrix<T>::FirstLocalRow() const
{ return graph_.FirstLocalSource(); }

template<typename T>
inline int
DistSparseMatrix<T>::LocalHeight() const
{ return graph_.NumLocalSources(); }

template<typename T>
inline int
DistSparseMatrix<T>::NumLocalEntries() const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::NumLocalEntries");
    EnsureConsistentSizes();
#endif
    return graph_.NumLocalEdges();
}

template<typename T>
inline int
DistSparseMatrix<T>::Capacity() const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
#endif
    return graph_.Capacity();
}

template<typename T>
inline int
DistSparseMatrix<T>::Row( int localIndex ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistSparseMatrix::Row");
#endif
    return graph_.Source( localIndex );
}

template<typename T>
inline int
DistSparseMatrix<T>::Col( int localIndex ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistSparseMatrix::Col");
#endif
    return graph_.Target( localIndex );
}

template<typename T>
inline int
DistSparseMatrix<T>::LocalEntryOffset( int localRow ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::LocalEntryOffset");
#endif
    return graph_.LocalEdgeOffset( localRow );
}

template<typename T>
inline int
DistSparseMatrix<T>::NumConnections( int localRow ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::NumConnections");
#endif
    return graph_.NumConnections( localRow );
}

template<typename T>
inline T
DistSparseMatrix<T>::Value( int localIndex ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistSparseMatrix::Value");
    if( localIndex < 0 || localIndex >= (int)values_.size() )
        throw std::logic_error("Entry number out of bounds");
#endif
    return values_[localIndex];
}

template<typename T>
inline bool
DistSparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{
    return a.i < b.i || (a.i == b.i && a.j < b.j);
}

template<typename T>
inline void
DistSparseMatrix<T>::StartAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::StartAssembly");
#endif
    graph_.EnsureNotAssembling();
    graph_.assembling_ = true;
}

template<typename T>
inline void
DistSparseMatrix<T>::StopAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::StopAssembly");
#endif
    if( !graph_.assembling_ )
        throw std::logic_error("Cannot stop assembly without starting");
    graph_.assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !graph_.sorted_ )
    {
        const int numLocalEntries = values_.size();
        std::vector<Entry<T> > entries( numLocalEntries );
        for( int s=0; s<numLocalEntries; ++s )
        {
            entries[s].i = graph_.sources_[s];
            entries[s].j = graph_.targets_[s];
            entries[s].value = values_[s];
        }
        std::sort( entries.begin(), entries.end(), CompareEntries );

        // Compress out duplicates
        int lastUnique=0;
        for( int s=1; s<numLocalEntries; ++s )
        {
            if( entries[s].i != entries[lastUnique].i ||
                entries[s].j != entries[lastUnique].j )
            {
                ++lastUnique;
                entries[lastUnique].i = entries[s].i;
                entries[lastUnique].j = entries[s].j;
                entries[lastUnique].value = entries[s].value;
            }
            else
                entries[lastUnique].value += entries[s].value;
        }
        const int numUnique = lastUnique+1;

        graph_.sources_.resize( numUnique );
        graph_.targets_.resize( numUnique );
        values_.resize( numUnique );
        for( int s=0; s<numUnique; ++s )
        {
            graph_.sources_[s] = entries[s].i;
            graph_.targets_[s] = entries[s].j;
            values_[s] = entries[s].value;
        }
    }
    graph_.ComputeLocalEdgeOffsets();
}

template<typename T>
inline void
DistSparseMatrix<T>::Reserve( int numLocalEntries )
{ 
    graph_.Reserve( numLocalEntries );
    values_.reserve( numLocalEntries );
}

template<typename T>
inline void
DistSparseMatrix<T>::Update( int row, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::Update");
    EnsureConsistentSizes();
#endif
    graph_.Insert( row, col );
    values_.push_back( value );
}

template<typename T>
inline void
DistSparseMatrix<T>::Empty()
{
    graph_.Empty();
    values_.clear();
}

template<typename T>
inline void
DistSparseMatrix<T>::ResizeTo( int height, int width )
{
    graph_.ResizeTo( height, width );
    values_.clear();
}

template<typename T>
inline void
DistSparseMatrix<T>::Print( std::string msg ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::Print");
#endif
    const int commSize = mpi::CommSize( graph_.comm_ );
    const int commRank = mpi::CommRank( graph_.comm_ );

    const int numLocalNonzeros = values_.size();
    std::vector<int> nonzeroSizes(commSize), nonzeroOffsets(commSize);
    mpi::AllGather( &numLocalNonzeros, 1, &nonzeroSizes[0], 1, graph_.comm_ );
    int numNonzeros=0;
    for( int q=0; q<commSize; ++q )
    {
        nonzeroOffsets[q] = numNonzeros;
        numNonzeros += nonzeroSizes[q];
    }

    std::vector<int> sources, targets;
    std::vector<T> values;
    if( commRank == 0 )
    {
        sources.resize( numNonzeros );
        targets.resize( numNonzeros );
        values.resize( numNonzeros );
    }
    mpi::Gather
    ( &graph_.sources_[0], numLocalNonzeros,
      &sources[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, graph_.comm_ );
    mpi::Gather
    ( &graph_.targets_[0], numLocalNonzeros,
      &targets[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, graph_.comm_ );
    mpi::Gather
    ( &values_[0], numLocalNonzeros,
      &values[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, graph_.comm_ );

    if( commRank == 0 )
    {
        if( msg != "" )
            std::cout << msg << std::endl;
        for( int s=0; s<numNonzeros; ++s )
            std::cout << sources[s] << " " 
                      << targets[s] << " " 
                      << values[s] << "\n";
        std::cout << std::endl;
    }
}

template<typename T>
inline void
DistSparseMatrix<T>::EnsureConsistentSizes() const
{ 
    graph_.EnsureConsistentSizes();
    if( graph_.NumLocalEdges() != (int)values_.size() )
        throw std::logic_error("Inconsistent sparsity sizes");
}

template<typename T>
inline void
DistSparseMatrix<T>::EnsureConsistentCapacities() const
{ 
    graph_.EnsureConsistentCapacities();
    if( graph_.Capacity() != values_.capacity() )
        throw std::logic_error("Inconsistent sparsity capacities");
}

} // namespace cliq
