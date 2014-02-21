/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
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
  SYMM_1D,                SYMM_2D,
  LDL_1D,                 LDL_2D,
  LDL_SELINV_1D,          LDL_SELINV_2D,
  LDL_INTRAPIV_1D,        LDL_INTRAPIV_2D,
  LDL_INTRAPIV_SELINV_1D, LDL_INTRAPIV_SELINV_2D,
  BLOCK_LDL_1D,           BLOCK_LDL_2D,
  BLOCK_LDL_INTRAPIV_1D,  BLOCK_LDL_INTRAPIV_2D
};

inline bool
Unfactored( SymmFrontType type )
{ return type == SYMM_1D || type == SYMM_2D; }

inline bool
FrontsAre1d( SymmFrontType type )
{
    return type == SYMM_1D                ||
           type == LDL_1D                 ||
           type == LDL_SELINV_1D          ||
           type == LDL_INTRAPIV_1D        ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == BLOCK_LDL_1D           ||
           type == BLOCK_LDL_INTRAPIV_1D;
}

inline bool
BlockFactorization( SymmFrontType type )
{ 
    return type == BLOCK_LDL_1D ||
           type == BLOCK_LDL_2D ||
           type == BLOCK_LDL_INTRAPIV_1D || 
           type == BLOCK_LDL_INTRAPIV_2D; 
}

inline bool
SelInvFactorization( SymmFrontType type )
{
    return type == LDL_SELINV_1D ||
           type == LDL_SELINV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D;
}

inline bool
PivotedFactorization( SymmFrontType type )
{
    return type == LDL_INTRAPIV_1D ||
           type == LDL_INTRAPIV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D ||
           type == BLOCK_LDL_INTRAPIV_1D  ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

inline SymmFrontType
ConvertTo2d( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo2d"))
    SymmFrontType newType;
    switch( type )
    {
    case SYMM_1D:                
    case SYMM_2D:                newType = SYMM_2D;                break;
    case LDL_1D: 
    case LDL_2D:                 newType = LDL_2D;                 break;
    case LDL_SELINV_1D:
    case LDL_SELINV_2D:          newType = LDL_SELINV_2D;          break;
    case LDL_INTRAPIV_1D:
    case LDL_INTRAPIV_2D:        newType = LDL_INTRAPIV_2D;        break;
    case LDL_INTRAPIV_SELINV_1D:
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_SELINV_2D; break;
    case BLOCK_LDL_1D:
    case BLOCK_LDL_2D:           newType = BLOCK_LDL_2D;           break;
    case BLOCK_LDL_INTRAPIV_1D:
    case BLOCK_LDL_INTRAPIV_2D:  newType = BLOCK_LDL_INTRAPIV_2D;  break;
    default: LogicError("Invalid front type");
    }
    return newType;
}

inline SymmFrontType
ConvertTo1d( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo1d"))
    SymmFrontType newType;
    switch( type )
    {
    case SYMM_1D:
    case SYMM_2D:                newType = SYMM_1D;                break;
    case LDL_1D:
    case LDL_2D:                 newType = LDL_1D;                 break;
    case LDL_SELINV_1D:
    case LDL_SELINV_2D:          newType = LDL_SELINV_1D;          break;
    case LDL_INTRAPIV_1D:
    case LDL_INTRAPIV_2D:        newType = LDL_INTRAPIV_1D;        break;
    case LDL_INTRAPIV_SELINV_1D:
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_SELINV_1D; break;
    case BLOCK_LDL_1D:
    case BLOCK_LDL_2D:           newType = BLOCK_LDL_1D;           break;
    case BLOCK_LDL_INTRAPIV_1D:
    case BLOCK_LDL_INTRAPIV_2D:  newType = BLOCK_LDL_INTRAPIV_1D;  break;
    default: LogicError("Invalid front type");
    }
    return newType;
}

inline SymmFrontType
AppendSelInv( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("AppendSelInv"))
    SymmFrontType newType;
    switch( type )
    {
    case LDL_1D:          newType = LDL_SELINV_1D; break;
    case LDL_2D:          newType = LDL_SELINV_2D; break;
    case LDL_INTRAPIV_1D: newType = LDL_INTRAPIV_SELINV_1D; break;
    case LDL_INTRAPIV_2D: newType = LDL_INTRAPIV_SELINV_2D; break;
    default: LogicError("Sel-inv does not make sense for this type");
    }
    return newType;
}

inline SymmFrontType
RemoveSelInv( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("RemoveSelInv"))
    SymmFrontType newType;
    switch( type )
    {
    case LDL_SELINV_1D: newType = LDL_1D; break;
    case LDL_SELINV_2D: newType = LDL_2D; break;
    case LDL_INTRAPIV_SELINV_1D: newType = LDL_INTRAPIV_1D; break;
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_2D; break;
    default: LogicError("This type did not involve selective inversion");
    }
    return newType;
}

inline SymmFrontType
InitialFactorType( SymmFrontType type )
{
    if( Unfactored(type) )
        LogicError("Front type does not require factorization");
    if( BlockFactorization(type) )
        return ConvertTo2d(type);
    else if( PivotedFactorization(type) )
        return LDL_INTRAPIV_2D;
    else
        return LDL_2D;
}

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename F>
struct SymmFront
{
    Matrix<F> frontL;

    Matrix<F> diag;
    Matrix<F> subdiag;
    Matrix<Int> piv;

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
    DistMatrix<F> front2dL;

    DistMatrix<F,VC,STAR> diag1d;
    DistMatrix<F,VC,STAR> subdiag1d;
    DistMatrix<Int,VC,STAR> piv;

    mutable DistMatrix<F,VC,STAR> work1d;
    mutable DistMatrix<F> work2d;
};

template<typename F>
struct DistSymmFrontTree
{
    bool isHermitian;
    SymmFrontType frontType;
    std::vector<SymmFront<F>> localFronts;
    std::vector<DistSymmFront<F>> distFronts;

    DistSymmFrontTree();

    DistSymmFrontTree
    ( const DistSparseMatrix<F>& A, 
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info,
      bool conjugate=false );

    void Initialize
    ( const DistSparseMatrix<F>& A,
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info,
      bool conjugate=false );

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
