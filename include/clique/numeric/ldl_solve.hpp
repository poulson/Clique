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
#ifndef CLIQUE_NUMERIC_LDL_SOLVE_HPP
#define CLIQUE_NUMERIC_LDL_SOLVE_HPP 1

namespace clique {
namespace numeric {

template<typename F>
void LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S, 
  const numeric::LocalSymmFact<F>& L, F alpha, Matrix<F>& X );

template<typename F>
void DistLDLForwardSolve
( const symbolic::LocalSymmFact& localS,
  const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<F>& localL,
  const numeric::DistSymmFact<F>& distL,
        Matrix<F>& localX );

template<typename F>
void LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<F>& L, Matrix<F>& X,
  bool checkIfSingular=false );

template<typename F>
void LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S, 
  const numeric::LocalSymmFact<F>& L, F alpha, Matrix<F>& X );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_LDL_SOLVE_HPP */

