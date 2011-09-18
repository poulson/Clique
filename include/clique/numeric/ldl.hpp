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
struct LocalSymmFact
{
    std::vector<elemental::Matrix<F> > fronts;
    mutable std::vector<elemental::Matrix<F> > solutions;
};

template<typename F>
struct DistSymmFact
{
    std::vector<elemental::DistMatrix<F,elemental::MC,elemental::MR> > fronts;
};

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling the following factorizations.

template<typename F>
void LocalLDL
( elemental::Orientation orientation, 
  symbolic::LocalSymmFact& S, // can't be const due to map...
  numeric::LocalSymmFact<F>& L );

template<typename F>
void DistLDL
( elemental::Orientation orientation,
        symbolic::DistSymmFact& S, // can't be const due to map...
  const numeric::LocalSymmFact<F>& localL,
        numeric::DistSymmFact<F>&  distL );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_LDL_HPP */

