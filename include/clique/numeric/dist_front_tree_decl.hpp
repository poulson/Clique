/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
struct Front
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
    std::vector<Front<F> > fronts;
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
