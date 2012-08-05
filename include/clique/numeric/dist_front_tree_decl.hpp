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

enum FrontType
{
  GENERAL_1D,
  GENERAL_2D,
  STRUCT_SYMM_1D,
  STRUCT_SYMM_2D,
  LU_1D,
  LU_2D,
  LU_SELINV_1D,
  LU_SELINV_2D,
  BLOCK_LU_2D
};

inline bool
FrontsAre1d( FrontType frontType )
{
    if( frontType == GENERAL_1D || 
        frontType == STRUCT_SYMM_1D || 
        frontType == LU_1D || 
        frontType == LU_SELINV_1D )
        return true;
    else
        return false;
}

template<typename F>
struct LocalFront
{
    Matrix<F> frontTL, frontBL, frontTR;
    mutable Matrix<F> work;
};

template<typename F>
struct DistFront
{
    // The 'frontType' member variable of the parent 'DistFrontTree' 
    // determines which of the following fronts is active.

    // TODO: Think about the fact that almost everything is now mutable...
    DistMatrix<int,VC,STAR> pivots;

    mutable DistMatrix<F,VC,STAR> front1dTL, front1dBL, front1dTR;
    mutable DistMatrix<F,VC,STAR> work1d;

    mutable DistMatrix<F> front2dTL, front2dBL, front2dTR;
    mutable DistMatrix<F> work2d;
};

template<typename F>
struct DistFrontTree
{
    FrontType frontType;
    std::vector<LocalFront<F> > localFronts;
    std::vector<DistFront<F> > distFronts;

    DistFrontTree();

    // For constructing structurally symmetric frontal matrices
    DistFrontTree
    ( const DistSparseMatrix<F>& A, 
      const DistMap& map,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info );
};

} // namespace cliq
