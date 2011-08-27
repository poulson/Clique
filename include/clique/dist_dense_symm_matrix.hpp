/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
#ifndef CLIQUE_DIST_DENSE_SYMM_MATRIX_HPP
#define CLIQUE_DIST_DENSE_SYMM_MATRIX_HPP 1

#include "mpi.h"

namespace clique {

// A simple symmetric 2d block-cyclic dense distributed matrix. Since it is for 
// internal usage only, we can require that the upper-left block is full and owned
// by the top-left process in the grid. We can also restrict access to blocks and 
// column panels of the lower triangle in order to facilitate packed storage, as
// well as require that the distribution blocks are square.
template<typename T>
class DistDenseSymmMatrix
{
private:
    int height_;
    int blockSize_;

    MPI_Comm comm_;
    int gridHeight_, gridWidth_;
    int gridRow_, gridCol_;

    std::vector<T> buffer_;

public:
    DistDenseSymmMatrix
    ( int height, int blockSize, MPI_Comm comm, int gridHeight, int gridWidth );
    ~DistDenseSymmMatrix();

    T* BlockBuffer( int i, int j );
    const T* BlockBuffer( int i, int j ) const;

    int BlockLDim( int i, int j );
};

} // namespace clique

#endif /* CLIQUE_DIST_DENSE_SYMM_MATRIX_HPP */
