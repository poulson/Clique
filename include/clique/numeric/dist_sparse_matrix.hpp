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
#ifndef CLIQUE_DIST_SPARSE_MATRIX_HPP
#define CLIQUE_DIST_SPARSE_MATRIX_HPP 1

namespace cliq {

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*ceil(height/commSize)
//     otherwise,        ceil(height/commSize)
template<typename F>
class DistSparseMatrix
{
public:
    DistSparseMatrix( int height, int width, mpi::Comm comm );

    int Height() const;
    int Width() const;
    const DistGraph& Graph() const;

    mpi::Comm Comm() const;
    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    int NumLocalEntries() const;
    int Capacity() const;

    void Reserve( int numLocalEntries );
    void PushBack( int i, int j, F value );

    void Empty();
    void ResizeTo( int height, int width );

private:
    DistGraph graph_;

    std::vector<F> entries_;

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline DistSparseMatrix<F>::DistSparseMatrix
( int height, int width, mpi::Comm comm )
: graph_(height,width,comm)
{ }

template<typename F>
inline int 
DistSparseMatrix<F>::Height() const
{ return graph_.NumSources(); }

template<typename F>
inline int 
DistSparseMatrix<F>::Width() const
{ return graph_.NumTargets(); }

template<typename F>
inline const DistGraph& 
DistSparseMatrix<F>::Graph() const
{ return graph_; }

template<typename F>
inline mpi::Comm 
DistSparseMatrix<F>::Comm() const
{ return graph_.Comm(); }

template<typename F>
inline int
DistSparseMatrix<F>::Blocksize() const
{ return graph_.Blocksize(); }

template<typename F>
inline int
DistSparseMatrix<F>::FirstLocalRow() const
{ return graph_.FirstLocalSource(); }

template<typename F>
inline int
DistSparseMatrix<F>::LocalHeight() const
{ return graph_.NumLocalSources(); }

template<typename F>
inline int
DistSparseMatrix<F>::NumLocalEntries() const
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::NumLocalEntries");
    EnsureConsistentSizes();
    PopCallStack();
#endif
    return graph_.NumLocalEdges();
}

template<typename F>
inline int
DistSparseMatrix<F>::Capacity() const
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
    PopCallStack();
#endif
    return graph_.Capacity();
}

template<typename F>
inline void
DistSparseMatrix<F>::Reserve( int numLocalEntries )
{ 
    graph_.Reserve( numLocalEntries );
    entries_.reserve( numLocalEntries );
}

template<typename F>
inline void
DistSparseMatrix<F>::PushBack( int i, int j, F value )
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::PushBack");
    EnsureConsistentSizes();
#endif
    graph_.PushBack( i, j );
    entries_.push_back( value );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistSparseMatrix<F>::Empty()
{
    graph_.Empty();
    entries_.clear();
}

template<typename F>
inline void
DistSparseMatrix<F>::ResizeTo( int height, int width )
{
    graph_.ResizeTo( height, width );
    entries_.clear();
}

template<typename F>
inline void
DistSparseMatrix<F>::EnsureConsistentSizes() const
{ 
    graph_.EnsureConsistentSizes();
    if( graph_.NumLocalEdges() != entries_.size() )
        throw std::logic_error("Inconsistent sparsity sizes");
}

template<typename F>
inline void
DistSparseMatrix<F>::EnsureConsistentCapacities() const
{ 
    graph_.EnsureConsistentCapacities();
    if( graph_.Capacity() != entries_.capacity() )
        throw std::logic_error("Inconsistent sparsity capacities");
}

} // namespace cliq

#endif // CLIQUE_DIST_SPARSE_MATRIX_HPP
