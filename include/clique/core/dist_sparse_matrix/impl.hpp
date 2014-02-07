/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_CORE_DISTSPARSEMATRIX_IMPL_HPP
#define CLIQ_CORE_DISTSPARSEMATRIX_IMPL_HPP

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
    SwapClear( vals_ );
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
    DEBUG_ONLY(
        CallStackEntry cse("DistSparseMatrix::NumLocalEntries");
        EnsureConsistentSizes();
    )
    return distGraph_.NumLocalEdges();
}

template<typename T>
inline int
DistSparseMatrix<T>::Capacity() const
{
    DEBUG_ONLY(
        CallStackEntry cse("DistSparseMatrix::Capacity");
        EnsureConsistentSizes();
        EnsureConsistentCapacities();
    )
    return distGraph_.Capacity();
}

template<typename T>
inline int
DistSparseMatrix<T>::Row( int localInd ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Row"))
    return distGraph_.Source( localInd );
}

template<typename T>
inline int
DistSparseMatrix<T>::Col( int localInd ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Col"))
    return distGraph_.Target( localInd );
}

template<typename T>
inline int
DistSparseMatrix<T>::LocalEntryOffset( int localRow ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::LocalEntryOffset"))
    return distGraph_.LocalEdgeOffset( localRow );
}

template<typename T>
inline int
DistSparseMatrix<T>::NumConnections( int localRow ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::NumConnections"))
    return distGraph_.NumConnections( localRow );
}

template<typename T>
inline T
DistSparseMatrix<T>::Value( int localInd ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("DistSparseMatrix::Value");
        if( localInd < 0 || localInd >= (int)vals_.size() )
            LogicError("Entry number out of bounds");
    )
    return vals_[localInd];
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
{ return &vals_[0]; }

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
{ return &vals_[0]; }

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
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::StartAssembly"))
    multMeta.ready = false;
    distGraph_.EnsureNotAssembling();
    distGraph_.assembling_ = true;
}

template<typename T>
inline void
DistSparseMatrix<T>::StopAssembly()
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::StopAssembly"))
    if( !distGraph_.assembling_ )
        LogicError("Cannot stop assembly without starting");
    distGraph_.assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !distGraph_.sorted_ )
    {
        const int numLocalEntries = vals_.size();
        std::vector<Entry<T> > entries( numLocalEntries );
        for( int s=0; s<numLocalEntries; ++s )
        {
            entries[s].i = distGraph_.sources_[s];
            entries[s].j = distGraph_.targets_[s];
            entries[s].value = vals_[s];
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
        vals_.resize( numUnique );
        for( int s=0; s<numUnique; ++s )
        {
            distGraph_.sources_[s] = entries[s].i;
            distGraph_.targets_[s] = entries[s].j;
            vals_[s] = entries[s].value;
        }
    }
    distGraph_.ComputeLocalEdgeOffsets();
}

template<typename T>
inline void
DistSparseMatrix<T>::Reserve( int numLocalEntries )
{ 
    distGraph_.Reserve( numLocalEntries );
    vals_.reserve( numLocalEntries );
}

template<typename T>
inline void
DistSparseMatrix<T>::Update( int row, int col, T value )
{
    DEBUG_ONLY(
        CallStackEntry cse("DistSparseMatrix::Update");
        EnsureConsistentSizes();
    )
    distGraph_.Insert( row, col );
    vals_.push_back( value );
}

template<typename T>
inline void
DistSparseMatrix<T>::Empty()
{
    distGraph_.Empty();
    SwapClear( vals_ );
}

template<typename T>
inline void
DistSparseMatrix<T>::Resize( int height, int width )
{
    distGraph_.Resize( height, width );
    SwapClear( vals_ );
}

template<typename T>
inline void
DistSparseMatrix<T>::EnsureConsistentSizes() const
{ 
    distGraph_.EnsureConsistentSizes();
    if( distGraph_.NumLocalEdges() != (int)vals_.size() )
        LogicError("Inconsistent sparsity sizes");
}

template<typename T>
inline void
DistSparseMatrix<T>::EnsureConsistentCapacities() const
{ 
    distGraph_.EnsureConsistentCapacities();
    if( distGraph_.Capacity() != vals_.capacity() )
        LogicError("Inconsistent sparsity capacities");
}

} // namespace cliq

#endif // ifndef CLIQ_CORE_DISTSPARSEMATRIX_IMPL_HPP
