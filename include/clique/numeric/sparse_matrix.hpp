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
#ifndef CLIQUE_SPARSE_MATRIX_HPP
#define CLIQUE_SPARSE_MATRIX_HPP 1

namespace cliq {

template<typename F>
class SparseMatrix
{
public:
    SparseMatrix();
    SparseMatrix( int height, int width );
    ~SparseMatrix();

    int Height() const;
    int Width() const;
    const cliq::Graph& Graph() const;

    int NumEntries() const;
    int Capacity() const;

    int Row( int entry ) const;
    int Col( int entry ) const;
    F Value( int entry ) const;
    int EntryOffset( int row ) const;
    int NumConnections( int row ) const;

    void Reserve( int numEntries );
    void PushBack( int row, int col, F value );

    void Empty();
    void ResizeTo( int height, int width );

private:
    cliq::Graph graph_;
    std::vector<F> values_;

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline
SparseMatrix<F>::SparseMatrix()
{ }

template<typename F>
inline 
SparseMatrix<F>::SparseMatrix( int height, int width )
: graph_(height,width)
{ }

template<typename F>
inline
SparseMatrix<F>::~SparseMatrix()
{ }

template<typename F>
inline int 
SparseMatrix<F>::Height() const
{ return graph_.NumSources(); }

template<typename F>
inline int 
SparseMatrix<F>::Width() const
{ return graph_.NumTargets(); }

template<typename F>
inline const cliq::Graph& 
SparseMatrix<F>::Graph() const
{ return graph_; }

template<typename F>
inline int
SparseMatrix<F>::NumEntries() const
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::NumEntries");
    EnsureConsistentSizes();
    PopCallStack();
#endif
    return graph_.NumEdges();
}

template<typename F>
inline int
SparseMatrix<F>::Capacity() const
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
    PopCallStack();
#endif
    return graph_.Capacity();
}

template<typename F>
inline int
SparseMatrix<F>::Row( int entry ) const
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

template<typename F>
inline int
SparseMatrix<F>::Col( int entry ) const
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

template<typename F>
inline int
SparseMatrix<F>::EntryOffset( int row ) const
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

template<typename F>
inline int
SparseMatrix<F>::NumConnections( int row ) const
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

template<typename F>
inline F
SparseMatrix<F>::Value( int entry ) const
{ 
#ifndef RELEASE 
    PushCallStack("SparseMatrix::Value");
    if( entry < 0 || entry >= values_.size() )
        throw std::logic_error("Entry number out of bounds");
    PopCallStack();
#endif
    return values_[entry];
}

template<typename F>
inline void
SparseMatrix<F>::Reserve( int numEntries )
{ 
    graph_.Reserve( numEntries );
    values_.reserve( numEntries );
}

template<typename F>
inline void
SparseMatrix<F>::PushBack( int row, int col, F value )
{
#ifndef RELEASE
    PushCallStack("SparseMatrix::PushBack");
    EnsureConsistentSizes();
#endif
    graph_.PushBack( row, col );
    values_.push_back( value );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SparseMatrix<F>::Empty()
{
    graph_.Empty();
    values_.clear();
}

template<typename F>
inline void
SparseMatrix<F>::ResizeTo( int height, int width )
{
    graph_.ResizeTo( height, width );
    values_.clear();
}

template<typename F>
inline void
SparseMatrix<F>::EnsureConsistentSizes() const
{ 
    graph_.EnsureConsistentSizes();
    if( graph_.NumEdges() != values_.size() )
        throw std::logic_error("Inconsistent sparsity sizes");
}

template<typename F>
inline void
SparseMatrix<F>::EnsureConsistentCapacities() const
{ 
    graph_.EnsureConsistentCapacities();
    if( graph_.Capacity() != values_.capacity() )
        throw std::logic_error("Inconsistent sparsity capacities");
}

} // namespace cliq

#endif // CLIQUE_SPARSE_MATRIX_HPP
