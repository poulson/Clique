/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename T>
inline
SparseMatrix<T>::SparseMatrix()
{ }

template<typename T>
inline
SparseMatrix<T>::SparseMatrix( int height )
: graph_(height)
{ }

template<typename T>
inline 
SparseMatrix<T>::SparseMatrix( int height, int width )
: graph_(height,width)
{ }

template<typename T>
inline
SparseMatrix<T>::SparseMatrix( const SparseMatrix<T>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::SparseMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct sparse matrix with itself");
}

template<typename T>
inline
SparseMatrix<T>::SparseMatrix( const DistSparseMatrix<T>& A )
{ 
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::SparseMatrix");
#endif
    *this = A;
}

template<typename T>
inline
SparseMatrix<T>::~SparseMatrix()
{ }

template<typename T>
inline int 
SparseMatrix<T>::Height() const
{ return graph_.NumSources(); }

template<typename T>
inline int 
SparseMatrix<T>::Width() const
{ return graph_.NumTargets(); }

template<typename T>
inline cliq::Graph& 
SparseMatrix<T>::Graph()
{ return graph_; }

template<typename T>
inline const cliq::Graph& 
SparseMatrix<T>::LockedGraph() const
{ return graph_; }

template<typename T>
inline int
SparseMatrix<T>::NumEntries() const
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::NumEntries");
    EnsureConsistentSizes();
#endif
    return graph_.NumEdges();
}

template<typename T>
inline int
SparseMatrix<T>::Capacity() const
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
#endif
    return graph_.Capacity();
}

template<typename T>
inline int
SparseMatrix<T>::Row( int index ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("SparseMatrix::Row");
#endif
    return graph_.Source( index );
}

template<typename T>
inline int
SparseMatrix<T>::Col( int index ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("SparseMatrix::Col");
#endif
    return graph_.Target( index );
}

template<typename T>
inline int
SparseMatrix<T>::EntryOffset( int row ) const
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::EntryOffset");
#endif
    return graph_.EdgeOffset( row );
}

template<typename T>
inline int
SparseMatrix<T>::NumConnections( int row ) const
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::NumConnections");
#endif
    return graph_.NumConnections( row );
}

template<typename T>
inline T
SparseMatrix<T>::Value( int index ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("SparseMatrix::Value");
    if( index < 0 || index >= values_.size() )
        throw std::logic_error("Entry number out of bounds");
#endif
    return values_[index];
}

template<typename T>
inline int*
SparseMatrix<T>::SourceBuffer()
{ return graph_.SourceBuffer(); }

template<typename T>
inline int*
SparseMatrix<T>::TargetBuffer()
{ return graph_.TargetBuffer(); }

template<typename T>
inline T*
SparseMatrix<T>::ValueBuffer()
{ return &values_[0]; }

template<typename T>
inline const int*
SparseMatrix<T>::LockedSourceBuffer() const
{ return graph_.LockedSourceBuffer(); }

template<typename T>
inline const int*
SparseMatrix<T>::LockedTargetBuffer() const
{ return graph_.LockedTargetBuffer(); }

template<typename T>
inline const T*
SparseMatrix<T>::LockedValueBuffer() const
{ return &values_[0]; }

template<typename T>
inline bool
SparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{ return a.i < b.i || (a.i == b.i && a.j < b.j); }

template<typename T>
inline void
SparseMatrix<T>::StartAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::StartAssembly");
#endif
    graph_.EnsureNotAssembling();
    graph_.assembling_ = true;
}

template<typename T>
inline void
SparseMatrix<T>::StopAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::StopAssembly");
#endif
    if( !graph_.assembling_ )
        throw std::logic_error("Cannot stop assembly without starting");
    graph_.assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !graph_.sorted_ )
    {
        const int numEntries = values_.size();
        std::vector<Entry<T> > entries( numEntries );
        for( int s=0; s<numEntries; ++s )
        {
            entries[s].i = graph_.sources_[s];
            entries[s].j = graph_.targets_[s];
            entries[s].value = values_[s];
        }
        std::sort( entries.begin(), entries.end(), CompareEntries );

        // Compress out duplicates
        int lastUnique=0;
        for( int s=1; s<numEntries; ++s )
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
    graph_.ComputeEdgeOffsets();
}

template<typename T>
inline void
SparseMatrix<T>::Reserve( int numEntries )
{ 
    graph_.Reserve( numEntries );
    values_.reserve( numEntries );
}

template<typename T>
inline void
SparseMatrix<T>::Update( int row, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::Update");
    EnsureConsistentSizes();
#endif
    graph_.Insert( row, col );
    values_.push_back( value );
}

template<typename T>
inline void
SparseMatrix<T>::Empty()
{
    graph_.Empty();
    values_.clear();
}

template<typename T>
inline void
SparseMatrix<T>::ResizeTo( int height, int width )
{
    graph_.ResizeTo( height, width );
    values_.clear();
}

template<typename T>
inline const SparseMatrix<T>&
SparseMatrix<T>::operator=( const SparseMatrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::operator=");
#endif
    graph_ = A.graph_;
    values_ = A.values_;
    return *this;
}

template<typename T>
inline const SparseMatrix<T>&
SparseMatrix<T>::operator=( const DistSparseMatrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SparseMatrix::operator=");
#endif
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    if( commSize != 1 )
        throw std::logic_error
        ("Can not yet construct from distributed sparse matrix");

    graph_ = A.graph_;
    values_ = A.values_;
    return *this;
}

template<typename T>
inline void
SparseMatrix<T>::EnsureConsistentSizes() const
{ 
    graph_.EnsureConsistentSizes();
    if( graph_.NumEdges() != values_.size() )
        throw std::logic_error("Inconsistent sparsity sizes");
}

template<typename T>
inline void
SparseMatrix<T>::EnsureConsistentCapacities() const
{ 
    graph_.EnsureConsistentCapacities();
    if( graph_.Capacity() != values_.capacity() )
        throw std::logic_error("Inconsistent sparsity capacities");
}

} // namespace cliq
