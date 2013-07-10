/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename F>
void Solve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F>
void Solve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void Solve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("Solve");
#endif
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D );
    const Orientation orientation = ( L.isHermitian ? ADJOINT : TRANSPOSE );
    if( !blockLDL )
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, L, X );
        // Solve against diagonal
        DiagonalSolve( info, L, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
    else
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, L, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
}

template<typename F>
inline void Solve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("Solve");
#endif
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D );
    const Orientation orientation = ( L.isHermitian ? ADJOINT : TRANSPOSE );
    if( !blockLDL )
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, info, L, X );
        // Solve against diagonal
        DiagonalSolve( info, L, X );
        // Solve against the (conjugate-)transpose of the unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
    else
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, info, L, X );
        // Solve against the (conjugate-)transpose of the block unit diagonal L
        LowerSolve( orientation, info, L, X );
    }
}

template<typename F>
inline void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, 
  bool sequential, int numDistSeps, int numSeqSeps, int cutoff )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricSolve");
#endif
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection
    ( A.LockedDistGraph(), map, sepTree, info, 
      sequential, numDistSeps, numSeqSeps, cutoff );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( TRANSPOSE, A, map, sepTree, info );
    LDL( info, frontTree, LDL_1D );

    DistNodalMultiVec<F> XNodal;
    XNodal.Pull( inverseMap, info, X );
    Solve( info, frontTree, XNodal );
    XNodal.Push( inverseMap, info, X );
}

template<typename F>
inline void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, 
  bool sequential, int numDistSeps, int numSeqSeps, int cutoff )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSolve");
#endif
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection
    ( A.LockedDistGraph(), map, sepTree, info, 
      sequential, numDistSeps, numSeqSeps, cutoff );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( ADJOINT, A, map, sepTree, info );
    LDL( info, frontTree, LDL_1D );

    DistNodalMultiVec<F> XNodal;
    XNodal.Pull( inverseMap, info, X );
    Solve( info, frontTree, XNodal );
    XNodal.Push( inverseMap, info, X );
}

} // namespace cliq
