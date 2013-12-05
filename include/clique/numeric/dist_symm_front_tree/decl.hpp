/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_DISTSYMMFRONTTREE_DECL_HPP
#define CLIQ_NUMERIC_DISTSYMMFRONTTREE_DECL_HPP

namespace cliq {

enum SymmFrontType
{
  SYMM_1D,
  SYMM_2D,
  LDL_1D,
  LDL_2D,
  LDL_SELINV_1D,
  LDL_SELINV_2D,
  BLOCK_LDL_2D
};

inline bool
FrontsAre1d( SymmFrontType frontType )
{
    return frontType == SYMM_1D ||
           frontType == LDL_1D  ||
           frontType == LDL_SELINV_1D;
}

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename F>
struct SymmFront
{
    Matrix<F> frontL;
    Matrix<F> diag;
    mutable Matrix<F> work;
};

template<typename F>
struct DistSymmFront
{
    // The 'frontType' member variable of the parent 'DistSymmFrontTree' 
    // determines which of the following fronts is active.
    //
    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    DistMatrix<F,VC,STAR> front1dL;
    DistMatrix<F,VC,STAR> diag1d;
    mutable DistMatrix<F,VC,STAR> work1d;

    DistMatrix<F> front2dL;
    mutable DistMatrix<F> work2d;
};

template<typename F>
struct DistSymmFrontTree
{
    bool isHermitian;
    SymmFrontType frontType;
    std::vector<SymmFront<F> > localFronts;
    std::vector<DistSymmFront<F> > distFronts;

    DistSymmFrontTree();

    DistSymmFrontTree
    ( Orientation orientation,
      const DistSparseMatrix<F>& A, 
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info );

    void Initialize
    ( Orientation orientation,
      const DistSparseMatrix<F>& A,
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info );

    void TopLeftMemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void BottomLeftMemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void MemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void FactorizationWork
    ( double& numLocalFlops, double& minLocalFlops, double& maxLocalFlops,
      double& numGlobalFlops, bool selInv=false ) const;

    void SolveWork
    ( double& numLocalFlops, double& minLocalFlops, double& maxLocalFlops,
      double& numGlobalFlops, int numRhs=1 ) const;
};

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_DISTSYMMFRONTTREE_DECL_HPP
