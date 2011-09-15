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
#ifndef CLIQUE_NUMERIC_SUPERNODE_LDL_HPP
#define CLIQUE_NUMERIC_SUPERNODE_LDL_HPP 1

namespace clique {
namespace numeric {

template<typename F>
void LocalSupernodeLDLT( elemental::Matrix<F>& A, int supernodeSize );

template<typename F>
void LocalSupernodeLDLH( elemental::Matrix<F>& A, int supernodeSize );

template<typename F>
void LocalSupernodeLDL
( elemental::Orientation orientation, 
  elemental::Matrix<F>& A, int supernodeSize );

template<typename F>
void DistSupernodeLDLT
( elemental::DistMatrix<F,elemental::MC,elemental::MR>& A, int supernodeSize );

template<typename F>
void DistSupernodeLDLH
( elemental::DistMatrix<F,elemental::MC,elemental::MR>& A, int supernodeSize );

template<typename F>
void DistSupernodeLDL
( elemental::Orientation orientation, 
  elemental::DistMatrix<F,elemental::MC,elemental::MR>& A, int supernodeSize );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_SUPERNODE_LDL_HPP */

