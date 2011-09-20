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

enum SolveMode { FEW_RHS, MANY_RHS };

namespace numeric {
using namespace elemental;

template<typename F>
struct LocalSymmFactSupernode
{
    Matrix<F> front;
    mutable Matrix<F> work;
};

template<typename F>
struct LocalSymmFact
{
    std::vector<LocalSymmFactSupernode<F> > supernodes;
};

template<typename F>
struct DistSymmFactSupernode
{
    // The 'SolveMode' member variable of the parent 'DistSymmFact' determines
    // which of the following fronts is active.
    //   FEW_RHS  -> front1d
    //   MANY_RHS -> front2d

    DistMatrix<F,VC,STAR> front1d;
    mutable DistMatrix<F,VC,STAR> work1d;

    DistMatrix<F,MC,MR> front2d;
    mutable DistMatrix<F,MC,MR> work2d;
};

template<typename F>
struct DistSymmFact
{
    SolveMode mode;
    std::vector<DistSymmFactSupernode<F> > supernodes;
};

template<typename F>
void SetSolveMode( DistSymmFact<F>& distL, SolveMode solveMode );

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling the following factorizations.

template<typename F>
void LocalLDL
( Orientation orientation, 
  symbolic::LocalSymmFact& S, // can't be const due to map...
  numeric::LocalSymmFact<F>& L );

template<typename F>
void DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S, // can't be const due to map...
  const numeric::LocalSymmFact<F>& localL,
        numeric::DistSymmFact<F>&  distL );

} // namespace numeric
} // namespace clique

#endif /* CLIQUE_NUMERIC_LDL_HPP */

