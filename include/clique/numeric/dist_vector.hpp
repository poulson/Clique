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
#ifndef CLIQUE_DIST_VECTOR_HPP
#define CLIQUE_DIST_VECTOR_HPP 1

namespace cliq {

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename F>
class DistVector
{
public:
    DistVector();
    DistVector( mpi::Comm comm );
    DistVector( int height, mpi::Comm comm );
    // TODO: Constructor for building from a DistVector
    ~DistVector();

    int Height() const;

    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    F GetLocal( int localRow ) const;
    void SetLocal( int localRow, F value );
    void UpdateLocal( int localRow, F value );

    void Empty();
    void ResizeTo( int height );

    const DistVector<F>& operator=( const DistVector<F>& x );

private:
    int height_;

    mpi::Comm comm_;

    int blocksize_;
    int firstLocalRow_;

    std::vector<F> values_;
};

// Set all of the entries of x to zero
template<typename F>
void MakeZeros( DistVector<F>& x );

// Draw the entries of x uniformly from the unitball in F
template<typename F>
void MakeUniform( DistVector<F>& x );

} // namespace cliq

#endif // CLIQUE_DIST_VECTOR_HPP
