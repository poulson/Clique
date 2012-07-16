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
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename F>
class DistSparseMatrix
{
public:
    DistSparseMatrix();
    DistSparseMatrix( mpi::Comm comm );
    DistSparseMatrix( int height, int width, mpi::Comm comm );
    // TODO: Constructor for building from another DistSparseMatrix
    ~DistSparseMatrix();

    int Height() const;
    int Width() const;
    const DistGraph& Graph() const;

    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    int NumLocalEntries() const;
    int Capacity() const;

    int Row( int localEntry ) const;
    int Col( int localEntry ) const;
    F Value( int localEntry ) const;
    int LocalEntryOffset( int localRow ) const;
    int NumConnections( int localRow ) const;

    void StartAssembly();
    void StopAssembly();
    void Reserve( int numLocalEntries );
    void PushBack( int row, int col, F value );

    void Empty();
    void ResizeTo( int height, int width );

    // TODO: operator=

private:
    DistGraph graph_;
    std::vector<F> values_;

    template<typename U>
    struct Entry
    {
        int i, j;
        U value;
    };

    static bool CompareEntries( const Entry<F>& a, const Entry<F>& b );

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    template<typename U> friend class SparseMatrix;
    template<typename U> friend class DistSymmFrontTree;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline 
DistSparseMatrix<F>::DistSparseMatrix()
{ }

template<typename F>
inline 
DistSparseMatrix<F>::DistSparseMatrix( mpi::Comm comm )
: graph_(comm)
{ }

template<typename F>
inline 
DistSparseMatrix<F>::DistSparseMatrix( int height, int width, mpi::Comm comm )
: graph_(height,width,comm)
{ }

template<typename F>
inline 
DistSparseMatrix<F>::~DistSparseMatrix()
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
inline void
DistSparseMatrix<F>::SetComm( mpi::Comm comm )
{ 
    graph_.SetComm( comm ); 
    values_.clear();
}

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
inline int
DistSparseMatrix<F>::Row( int localEntry ) const
{ 
#ifndef RELEASE 
    PushCallStack("DistSparseMatrix::Row");
#endif
    const int row = graph_.Source( localEntry );
#ifndef RELEASE
    PopCallStack();
#endif
    return row;
}

template<typename F>
inline int
DistSparseMatrix<F>::Col( int localEntry ) const
{ 
#ifndef RELEASE 
    PushCallStack("DistSparseMatrix::Col");
#endif
    const int col = graph_.Target( localEntry );
#ifndef RELEASE
    PopCallStack();
#endif
    return col;
}

template<typename F>
inline int
DistSparseMatrix<F>::LocalEntryOffset( int localRow ) const
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::LocalEntryOffset");
#endif
    const int localEntryOffset = graph_.LocalEdgeOffset( localRow );
#ifndef RELEASE
    PopCallStack();
#endif
    return localEntryOffset;
}

template<typename F>
inline int
DistSparseMatrix<F>::NumConnections( int localRow ) const
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::NumConnections");
#endif
    const int numConnections = graph_.NumConnections( localRow );
#ifndef RELEASE
    PopCallStack();
#endif
    return numConnections;
}

template<typename F>
inline F
DistSparseMatrix<F>::Value( int localEntry ) const
{ 
#ifndef RELEASE 
    PushCallStack("DistSparseMatrix::Value");
    if( localEntry < 0 || localEntry >= (int)values_.size() )
        throw std::logic_error("Entry number out of bounds");
    PopCallStack();
#endif
    return values_[localEntry];
}

template<typename F>
inline bool
DistSparseMatrix<F>::CompareEntries( const Entry<F>& a, const Entry<F>& b )
{
    return a.i < b.i || (a.i == b.i && a.j < b.j);
}

template<typename F>
inline void
DistSparseMatrix<F>::StartAssembly()
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::StartAssembly");
#endif
    graph_.EnsureNotAssembling();
    graph_.assembling_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistSparseMatrix<F>::StopAssembly()
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::StopAssembly");
#endif
    if( !graph_.assembling_ )
        throw std::logic_error("Cannot stop assembly without starting");
    graph_.assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !graph_.sorted_ )
    {
        const int numLocalEntries = values_.size();
        std::vector<Entry<F> > entries( numLocalEntries );
        for( int s=0; s<numLocalEntries; ++s )
        {
            entries[s].i = graph_.sources_[s];
            entries[s].j = graph_.targets_[s];
            entries[s].value = values_[s];
        }
        std::sort( entries.begin(), entries.end(), CompareEntries );
        for( int s=0; s<numLocalEntries; ++s )
        {
            graph_.sources_[s] = entries[s].i;
            graph_.targets_[s] = entries[s].j;
            values_[s] = entries[s].value;
        }
    }

    graph_.ComputeLocalEdgeOffsets();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistSparseMatrix<F>::Reserve( int numLocalEntries )
{ 
    graph_.Reserve( numLocalEntries );
    values_.reserve( numLocalEntries );
}

template<typename F>
inline void
DistSparseMatrix<F>::PushBack( int row, int col, F value )
{
#ifndef RELEASE
    PushCallStack("DistSparseMatrix::PushBack");
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
DistSparseMatrix<F>::Empty()
{
    graph_.Empty();
    values_.clear();
}

template<typename F>
inline void
DistSparseMatrix<F>::ResizeTo( int height, int width )
{
    graph_.ResizeTo( height, width );
    values_.clear();
}

template<typename F>
inline void
DistSparseMatrix<F>::EnsureConsistentSizes() const
{ 
    graph_.EnsureConsistentSizes();
    if( graph_.NumLocalEdges() != (int)values_.size() )
        throw std::logic_error("Inconsistent sparsity sizes");
}

template<typename F>
inline void
DistSparseMatrix<F>::EnsureConsistentCapacities() const
{ 
    graph_.EnsureConsistentCapacities();
    if( graph_.Capacity() != values_.capacity() )
        throw std::logic_error("Inconsistent sparsity capacities");
}

} // namespace cliq

#endif // CLIQUE_DIST_SPARSE_MATRIX_HPP
