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
#ifndef CLIQUE_NUMERIC_SUPERNODE_SOLVE_HPP
#define CLIQUE_NUMERIC_SUPERNODE_SOLVE_HPP 1

namespace clique {
namespace numeric {
using namespace elemental;

template<typename F>
void LocalSupernodeLDLForwardSolve
( int supernodeSize,
  const Matrix<F>& L, Matrix<F>& X );

template<typename F>
void LocalSupernodeLDLDiagonalSolve
( int supernodeSize,
  const Matrix<F>& d, Matrix<F>& X,
  bool checkIfSingular=false );

template<typename F>
void LocalSupernodeLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const Matrix<F>& U, Matrix<F>& X );

template<typename F>
void DistSupernodeLDLForwardSolve
( int supernodeSize,
  const DistMatrix<F,VC,STAR>& L, 
        DistMatrix<F,VC,STAR>& X );

template<typename F>
void DistSupernodeLDLDiagonalSolve
( int supernodeSize,
  const DistMatrix<F,VC,STAR>& d,
        DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular=false );

template<typename F>
void DistSupernodeLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const DistMatrix<F,VC,STAR>& U, 
        DistMatrix<F,VC,STAR>& X );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_SUPERNODE_SOLVE_HPP */

