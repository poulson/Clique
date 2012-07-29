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
    const DistGraph& Graph() const;

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
    int Row( int localEntry ) const;
    int Col( int localEntry ) const;
    T Value( int localEntry ) const;
    int NumLocalEntries() const;
    int LocalEntryOffset( int localRow ) const;
    int NumConnections( int localRow ) const;

    // For modifying the size of the matrix
    void Empty();
    void ResizeTo( int height, int width );

    // TODO: operator=

private:
    DistGraph graph_;
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
    template<typename U> friend class DistSymmFrontTree;
};


} // namespace cliq
