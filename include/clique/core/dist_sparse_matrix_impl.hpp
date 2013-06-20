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
DistSparseMatrix<T>::DistSparseMatrix()
{ }

template<typename T>
inline 
DistSparseMatrix<T>::DistSparseMatrix( mpi::Comm comm )
: distGraph_(comm)
{ }

template<typename T>
inline
DistSparseMatrix<T>::DistSparseMatrix( int height, mpi::Comm comm )
: distGraph_(height,comm)
{ }

template<typename T>
inline 
DistSparseMatrix<T>::DistSparseMatrix( int height, int width, mpi::Comm comm )
: distGraph_(height,width,comm)
{ }

template<typename T>
inline 
DistSparseMatrix<T>::~DistSparseMatrix()
{ }

template<typename T>
inline int 
DistSparseMatrix<T>::Height() const
{ return distGraph_.NumSources(); }

template<typename T>
inline int 
DistSparseMatrix<T>::Width() const
{ return distGraph_.NumTargets(); }

template<typename T>
inline cliq::DistGraph& 
DistSparseMatrix<T>::DistGraph()
{ return distGraph_; }

template<typename T>
inline const cliq::DistGraph& 
DistSparseMatrix<T>::LockedDistGraph() const
{ return distGraph_; }

template<typename T>
inline void
DistSparseMatrix<T>::SetComm( mpi::Comm comm )
{ 
    distGraph_.SetComm( comm ); 
    values_.clear();
}

template<typename T>
inline mpi::Comm 
DistSparseMatrix<T>::Comm() const
{ return distGraph_.Comm(); }

template<typename T>
inline int
DistSparseMatrix<T>::Blocksize() const
{ return distGraph_.Blocksize(); }

template<typename T>
inline int
DistSparseMatrix<T>::FirstLocalRow() const
{ return distGraph_.FirstLocalSource(); }

template<typename T>
inline int
DistSparseMatrix<T>::LocalHeight() const
{ return distGraph_.NumLocalSources(); }

template<typename T>
inline int
DistSparseMatrix<T>::NumLocalEntries() const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::NumLocalEntries");
    EnsureConsistentSizes();
#endif
    return distGraph_.NumLocalEdges();
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
    return distGraph_.Capacity();
}

template<typename T>
inline int
DistSparseMatrix<T>::Row( int localIndex ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistSparseMatrix::Row");
#endif
    return distGraph_.Source( localIndex );
}

template<typename T>
inline int
DistSparseMatrix<T>::Col( int localIndex ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistSparseMatrix::Col");
#endif
    return distGraph_.Target( localIndex );
}

template<typename T>
inline int
DistSparseMatrix<T>::LocalEntryOffset( int localRow ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::LocalEntryOffset");
#endif
    return distGraph_.LocalEdgeOffset( localRow );
}

template<typename T>
inline int
DistSparseMatrix<T>::NumConnections( int localRow ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::NumConnections");
#endif
    return distGraph_.NumConnections( localRow );
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
inline int*
DistSparseMatrix<T>::SourceBuffer()
{ return distGraph_.SourceBuffer(); }

template<typename T>
inline int*
DistSparseMatrix<T>::TargetBuffer()
{ return distGraph_.TargetBuffer(); }

template<typename T>
inline T*
DistSparseMatrix<T>::ValueBuffer()
{ return &values_[0]; }

template<typename T>
inline const int*
DistSparseMatrix<T>::LockedSourceBuffer() const
{ return distGraph_.LockedSourceBuffer(); }

template<typename T>
inline const int*
DistSparseMatrix<T>::LockedTargetBuffer() const
{ return distGraph_.LockedTargetBuffer(); }

template<typename T>
inline const T*
DistSparseMatrix<T>::LockedValueBuffer() const
{ return &values_[0]; }

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
    distGraph_.EnsureNotAssembling();
    distGraph_.assembling_ = true;
}

template<typename T>
inline void
DistSparseMatrix<T>::StopAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("DistSparseMatrix::StopAssembly");
#endif
    if( !distGraph_.assembling_ )
        throw std::logic_error("Cannot stop assembly without starting");
    distGraph_.assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !distGraph_.sorted_ )
    {
        const int numLocalEntries = values_.size();
        std::vector<Entry<T> > entries( numLocalEntries );
        for( int s=0; s<numLocalEntries; ++s )
        {
            entries[s].i = distGraph_.sources_[s];
            entries[s].j = distGraph_.targets_[s];
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

        distGraph_.sources_.resize( numUnique );
        distGraph_.targets_.resize( numUnique );
        values_.resize( numUnique );
        for( int s=0; s<numUnique; ++s )
        {
            distGraph_.sources_[s] = entries[s].i;
            distGraph_.targets_[s] = entries[s].j;
            values_[s] = entries[s].value;
        }
    }
    distGraph_.ComputeLocalEdgeOffsets();
}

template<typename T>
inline void
DistSparseMatrix<T>::Reserve( int numLocalEntries )
{ 
    distGraph_.Reserve( numLocalEntries );
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
    distGraph_.Insert( row, col );
    values_.push_back( value );
}

template<typename T>
inline void
DistSparseMatrix<T>::Empty()
{
    distGraph_.Empty();
    values_.clear();
}

template<typename T>
inline void
DistSparseMatrix<T>::ResizeTo( int height, int width )
{
    distGraph_.ResizeTo( height, width );
    values_.clear();
}

template<typename T>
inline void
DistSparseMatrix<T>::EnsureConsistentSizes() const
{ 
    distGraph_.EnsureConsistentSizes();
    if( distGraph_.NumLocalEdges() != (int)values_.size() )
        throw std::logic_error("Inconsistent sparsity sizes");
}

template<typename T>
inline void
DistSparseMatrix<T>::EnsureConsistentCapacities() const
{ 
    distGraph_.EnsureConsistentCapacities();
    if( distGraph_.Capacity() != values_.capacity() )
        throw std::logic_error("Inconsistent sparsity capacities");
}

} // namespace cliq
