/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_CORE_DISTSPARSEMATRIX_DECL_HPP
#define CLIQ_CORE_DISTSPARSEMATRIX_DECL_HPP

namespace cliq {

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename T>
class DistSparseMatrix
{
public:
    // Construction and destruction
    DistSparseMatrix();
    DistSparseMatrix( mpi::Comm comm );
    DistSparseMatrix( int height, mpi::Comm comm );
    DistSparseMatrix( int height, int width, mpi::Comm comm );
    // TODO: Constructor for building from another DistSparseMatrix
    ~DistSparseMatrix();

    // High-level information
    int Height() const;
    int Width() const;
    cliq::DistGraph& DistGraph();
    const cliq::DistGraph& LockedDistGraph() const;

    // Communicator-management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    // Assembly-related routines
    void StartAssembly();
    void StopAssembly();
    void Reserve( int numLocalEntries );
    void Update( int row, int col, T value );
    int Capacity() const;

    // Local data
    int Row( int localIndex ) const;
    int Col( int localIndex ) const;
    T Value( int localIndex ) const;
    int NumLocalEntries() const;
    int LocalEntryOffset( int localRow ) const;
    int NumConnections( int localRow ) const;

    int* SourceBuffer();
    int* TargetBuffer();
    T* ValueBuffer();
    const int* LockedSourceBuffer() const;
    const int* LockedTargetBuffer() const;
    const T* LockedValueBuffer() const;

    // For modifying the size of the matrix
    void Empty();
    void ResizeTo( int height, int width );

    // TODO: operator=

private:
    cliq::DistGraph distGraph_;
    std::vector<T> values_;

    template<typename U>
    struct Entry
    {
        int i, j;
        U value;
    };

    static bool CompareEntries( const Entry<T>& a, const Entry<T>& b );

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    template<typename U> friend class SparseMatrix;
    template<typename U> friend struct DistSymmFrontTree;
};

} // namespace cliq

#endif // ifndef CLIQ_CORE_DISTSPARSEMATRIX_DECL_HPP
