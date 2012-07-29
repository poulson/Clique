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
class DistMultiVector
{
public:
    // Constructors and destructors
    DistMultiVector();
    DistMultiVector( mpi::Comm comm );
    DistMultiVector( int height, int width, mpi::Comm comm );
    // TODO: Constructor for building from a DistVector
    ~DistMultiVector();

    // High-level information
    int Height() const;
    int Width() const;

    // Communicator management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    // Local data
    T GetLocal( int localRow, int col ) const;
    void SetLocal( int localRow, int col, T value );
    void UpdateLocal( int localRow, int col, T value );

    // For modifying the size of the multi-vector
    void Empty();
    void ResizeTo( int height, int width );

    // Assignment
    const DistMultiVector<T>& operator=( const DistVector<T>& x );
    const DistMultiVector<T>& operator=( const DistMultiVector<T>& X );

private:
    int height_, width_;

    mpi::Comm comm_;

    int blocksize_;
    int firstLocalRow_;

    Matrix<T> localMultiVec_;
};

// Set all of the entries of X to zero
template<typename T>
void MakeZeros( DistMultiVector<T>& X );

// Draw the entries of X uniformly from the unitball in T
template<typename T>
void MakeUniform( DistMultiVector<T>& X );

// Just an l2 norm for now
template<typename F>
void Norms
( const DistMultiVector<F>& X, std::vector<typename Base<F>::type>& norms );

// Y := alpha X + Y
template<typename T>
void Axpy( T alpha, const DistMultiVector<T>& X, DistMultiVector<T>& Y );

} // namespace cliq
