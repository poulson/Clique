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

template<typename F>
void Solve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistVector<F>& x,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistVector<F>& x,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );

template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVector<F>& X,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVector<F>& X,
  bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void Solve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("Solve");
#endif
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D );
    if( !blockLDL )
    {
        // Solve against unit diagonal L
        LowerSolve( NORMAL, UNIT, info, L, localX );

        // Solve against diagonal
        DiagonalSolve( info, L, localX );

        // Solve against the (conjugate-)transpose of the unit diagonal L
        if( L.isHermitian )
            LowerSolve( ADJOINT, UNIT, info, L, localX );
        else
            LowerSolve( TRANSPOSE, UNIT, info, L, localX );
    }
    else
    {
        // Solve against block diagonal factor, L D
        LowerSolve( NORMAL, NON_UNIT, info, L, localX );

        // Solve against the (conjugate-)transpose of the block unit diagonal L
        if( L.isHermitian )
            LowerSolve( ADJOINT, NON_UNIT, info, L, localX );
        else
            LowerSolve( TRANSPOSE, NON_UNIT, info, L, localX );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void SymmetricSolve
( const DistSparseMatrix<F>& A, DistVector<F>& x, 
  bool sequential, int numDistSeps, int numSeqSeps, int cutoff )
{
#ifndef RELEASE
    PushCallStack("SymmetricSolve");
#endif
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection
    ( A.Graph(), map, sepTree, info, 
      sequential, numDistSeps, numSeqSeps, cutoff );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( TRANSPOSE, A, map, sepTree, info );
    LDL( info, frontTree, LDL_1D );

    DistNodalVector<F> xNodal;
    xNodal.Pull( inverseMap, info, x );
    Solve( info, frontTree, xNodal.localVec );
    xNodal.Push( inverseMap, info, x );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSolve
( const DistSparseMatrix<F>& A, DistVector<F>& x, 
  bool sequential, int numDistSeps, int numSeqSeps, int cutoff )
{
#ifndef RELEASE
    PushCallStack("HermitianSolve");
#endif
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection
    ( A.Graph(), map, sepTree, info, 
      sequential, numDistSeps, numSeqSeps, cutoff );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( ADJOINT, A, map, sepTree, info );
    LDL( info, frontTree, LDL_1D );

    DistNodalVector<F> xNodal;
    xNodal.Pull( inverseMap, info, x );
    Solve( info, frontTree, xNodal.localVec );
    xNodal.Push( inverseMap, info, x );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVector<F>& X, 
  bool sequential, int numDistSeps, int numSeqSeps, int cutoff )
{
#ifndef RELEASE
    PushCallStack("SymmetricSolve");
#endif
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection
    ( A.Graph(), map, sepTree, info, 
      sequential, numDistSeps, numSeqSeps, cutoff );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( TRANSPOSE, A, map, sepTree, info );
    LDL( info, frontTree, LDL_1D );

    DistNodalMultiVector<F> XNodal;
    XNodal.Pull( inverseMap, info, X );
    Solve( info, frontTree, XNodal.localMultiVec );
    XNodal.Push( inverseMap, info, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVector<F>& X, 
  bool sequential, int numDistSeps, int numSeqSeps, int cutoff )
{
#ifndef RELEASE
    PushCallStack("HermitianSolve");
#endif
    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, inverseMap;
    NestedDissection
    ( A.Graph(), map, sepTree, info, 
      sequential, numDistSeps, numSeqSeps, cutoff );
    map.FormInverse( inverseMap );

    DistSymmFrontTree<F> frontTree( ADJOINT, A, map, sepTree, info );
    LDL( info, frontTree, LDL_1D );

    DistNodalMultiVector<F> XNodal;
    XNodal.Pull( inverseMap, info, X );
    Solve( info, frontTree, XNodal.localMultiVec );
    XNodal.Push( inverseMap, info, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
