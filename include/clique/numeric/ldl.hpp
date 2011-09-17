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
#ifndef CLIQUE_NUMERIC_LDL_HPP
#define CLIQUE_NUMERIC_LDL_HPP 1

namespace clique {
namespace numeric {

template<typename F>
struct LocalFactMatrix
{
    std::vector<elemental::Matrix<F> > fronts;
    mutable std::vector<elemental::Matrix<F> > solutions;
};

template<typename F>
struct DistFactMatrix
{
    std::vector<elemental::DistMatrix<F,elemental::MC,elemental::MR> > fronts;
};

// When this is called, all of the fronts in L should by equal to the expansions
// of the original sparse matrix.
template<typename F>
void LocalLDL
( elemental::Orientation orientation, 
  symbolic::LocalFactStruct& S, // can't be const due to map...
  numeric::LocalFactMatrix<F>& L );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_LDL_HPP */

