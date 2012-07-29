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
class DistVector
{
public:
    // Constructors and destructors
    DistVector();
    DistVector( mpi::Comm comm );
    DistVector( int height, mpi::Comm comm );
    // TODO: Constructor for building from a DistVector
    ~DistVector();

    // High-level information
    int Height() const;

    // Communicator management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    // Local data
    T GetLocal( int localRow ) const;
    void SetLocal( int localRow, T value );
    void UpdateLocal( int localRow, T value );

    // For modifying the size of the vector
    void Empty();
    void ResizeTo( int height );

    // Assignment
    const DistVector<T>& operator=( const DistVector<T>& x );

private:
    int height_;

    mpi::Comm comm_;

    int blocksize_;
    int firstLocalRow_;

    Matrix<T> localVec_;
};

// Set all of the entries of x to zero
template<typename T>
void MakeZeros( DistVector<T>& x );

// Draw the entries of x uniformly from the unitball in T
template<typename T>
void MakeUniform( DistVector<T>& x );

// Just an l2 norm for now
template<typename F>
typename Base<F>::type Norm( const DistVector<F>& x );

// y := alpha x + y
template<typename T>
void Axpy( T alpha, const DistVector<T>& x, DistVector<T>& y );

} // namespace cliq
