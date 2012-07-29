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
    PushCallStack("SparseMatrix::SparseMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct sparse matrix with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
SparseMatrix<T>::SparseMatrix( const DistSparseMatrix<T>& A )
{ 
#ifndef RELEASE
    PushCallStack("SparseMatrix::SparseMatrix");
#endif
    *this = A;
#ifndef RELEASE
    PopCallStack();
#endif
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
inline const cliq::Graph& 
SparseMatrix<T>::Graph() const
{ return graph_; }

template<typename T>
inline int
SparseMatrix<T>::NumEntries() const
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::NumEntries");
    EnsureConsistentSizes();
    PopCallStack();
#endif
    return graph_.NumEdges();
}

template<typename T>
inline int
SparseMatrix<T>::Capacity() const
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
    PopCallStack();
#endif
    return graph_.Capacity();
}

template<typename T>
inline int
SparseMatrix<T>::Row( int entry ) const
{ 
#ifndef RELEASE 
    PushCallStack("SparseMatrix::Row");
#endif
    const int row = graph_.Source( entry );
#ifndef RELEASE
    PopCallStack();
#endif
    return row;
}

template<typename T>
inline int
SparseMatrix<T>::Col( int entry ) const
{ 
#ifndef RELEASE 
    PushCallStack("SparseMatrix::Col");
#endif
    const int col = graph_.Target( entry );
#ifndef RELEASE
    PopCallStack();
#endif
    return col;
}

template<typename T>
inline int
SparseMatrix<T>::EntryOffset( int row ) const
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::EntryOffset");
#endif
    const int entryOffset = graph_.EdgeOffset( row );
#ifndef RELEASE
    PopCallStack();
#endif
    return entryOffset;
}

template<typename T>
inline int
SparseMatrix<T>::NumConnections( int row ) const
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::NumConnections");
#endif
    const int numConnections = graph_.NumConnections( row );
#ifndef RELEASE
    PopCallStack();
#endif
    return numConnections;
}

template<typename T>
inline T
SparseMatrix<T>::Value( int entry ) const
{ 
#ifndef RELEASE 
    PushCallStack("SparseMatrix::Value");
    if( entry < 0 || entry >= values_.size() )
        throw std::logic_error("Entry number out of bounds");
    PopCallStack();
#endif
    return values_[entry];
}

template<typename T>
inline bool
SparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{
    return a.i < b.i || (a.i == b.i && a.j < b.j);
}

template<typename T>
inline void
SparseMatrix<T>::StartAssembly()
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::StartAssembly");
#endif
    graph_.EnsureNotAssembling();
    graph_.assembling_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SparseMatrix<T>::StopAssembly()
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::StopAssembly");
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
#ifndef RELEASE
    PopCallStack();
#endif
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
    PushCallStack("SparseMatrix::Update");
    EnsureConsistentSizes();
#endif
    graph_.Insert( row, col );
    values_.push_back( value );
#ifndef RELEASE
    PopCallStack();
#endif
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
    PushCallStack("SparseMatrix::operator=");
#endif
    graph_ = A.graph_;
    values_ = A.values_;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template<typename T>
inline const SparseMatrix<T>&
SparseMatrix<T>::operator=( const DistSparseMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::operator=");
#endif
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    if( commSize != 1 )
        throw std::logic_error
        ("Can not yet construct from distributed sparse matrix");

    graph_ = A.graph_;
    values_ = A.values_;
#ifndef RELEASE
    PopCallStack();
#endif
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
