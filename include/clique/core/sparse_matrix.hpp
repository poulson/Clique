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

template<typename T>
class SparseMatrix
{
public:
    // Construction and destruction
    SparseMatrix();
    SparseMatrix( int height );
    SparseMatrix( int height, int width );
    SparseMatrix( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    SparseMatrix( const DistSparseMatrix<T>& A );
    ~SparseMatrix();

    // High-level information
    int Height() const;
    int Width() const;
    const cliq::Graph& Graph() const;

    // Assembly-related routines
    void StartAssembly();
    void StopAssembly();
    void Reserve( int numEntries );
    void Update( int row, int col, T value );
    int Capacity() const;

    // Data
    int Row( int entry ) const;
    int Col( int entry ) const;
    T Value( int entry ) const;
    int NumEntries() const;
    int EntryOffset( int row ) const;
    int NumConnections( int row ) const;

    // For modifying the size of the matrix
    void Empty();
    void ResizeTo( int height, int width );

    // For copying one matrix to another
    const SparseMatrix<T>& operator=( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    const SparseMatrix<T>& operator=( const DistSparseMatrix<T>& A );

private:
    cliq::Graph graph_;
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

    template<typename U> friend class DistSparseMatrix;
};

} // namespace cliq

#endif // CLIQUE_SPARSE_MATRIX_HPP
