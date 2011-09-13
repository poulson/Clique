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

template<typename F>
void SupernodeLDLForwardSolve
( int supernodeSize,
  F alpha, const elemental::DistMatrix<F,elemental::VC,elemental::STAR>& L, 
                 elemental::DistMatrix<F,elemental::VC,elemental::STAR>& X );

template<typename F>
void SupernodeLDLDiagonalSolve
( int supernodeSize,
  F alpha, const elemental::DistMatrix<F,elemental::VC,elemental::STAR>& d,
                 elemental::DistMatrix<F,elemental::VC,elemental::STAR>& X,
  bool checkIfSingular=false );

template<typename F>
void SupernodeLDLBackwardSolve
( elemental::Orientation orientation, int supernodeSize,
  F alpha, const elemental::DistMatrix<F,elemental::VC,elemental::STAR>& U, 
                 elemental::DistMatrix<F,elemental::VC,elemental::STAR>& X );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_SUPERNODE_SOLVE_HPP */

