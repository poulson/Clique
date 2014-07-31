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
#ifndef CLIQ_NUMERIC_DIAGONALSOLVE_HPP
#define CLIQ_NUMERIC_DIAGONALSOLVE_HPP

namespace cliq {

template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X );
template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X );
template<typename F>
void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X );

template<typename F>
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X );
template<typename F>
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X );

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalDiagonalSolve"))
    const Int numLocalNodes = info.localNodes.size();
    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=0; s<numLocalNodes; ++s )
            El::QuasiDiagonalSolve
            ( LEFT, LOWER, L.localFronts[s].diag, L.localFronts[s].subdiag, 
              X.localNodes[s], L.isHermitian );
    }
    else
    {
        for( Int s=0; s<numLocalNodes; ++s )
            El::DiagonalSolve
            ( LEFT, NORMAL, L.localFronts[s].diag, X.localNodes[s], true );
    }
}

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalDiagonalSolve"))
    const Int numLocalNodes = info.localNodes.size();
    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=0; s<numLocalNodes; ++s ) 
            El::QuasiDiagonalSolve
            ( LEFT, LOWER, L.localFronts[s].diag, L.localFronts[s].subdiag, 
              X.localNodes[s], L.isHermitian );
    }
    else
    {
        for( Int s=0; s<numLocalNodes; ++s )
            El::DiagonalSolve
            ( LEFT, NORMAL, L.localFronts[s].diag, X.localNodes[s], true );
    }
}

template<typename F> 
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistDiagonalSolve"))
    const Int numDistNodes = info.distNodes.size();

    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            El::QuasiDiagonalSolve
            ( LEFT, LOWER, front.diag1d, front.subdiag1d, X.distNodes[s-1],
              L.isHermitian );
        }
    }
    else
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            El::DiagonalSolve
            ( LEFT, NORMAL, front.diag1d, X.distNodes[s-1], true );
        }
    }
}

template<typename F> 
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistDiagonalSolve"))
    const Int numDistNodes = info.distNodes.size();

    if( PivotedFactorization(L.frontType) )
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            El::QuasiDiagonalSolve
            ( LEFT, LOWER, front.diag1d, front.subdiag1d, X.distNodes[s-1],
              L.isHermitian );
        }
    }
    else
    {
        for( Int s=1; s<numDistNodes; ++s )
        {
            const DistSymmFront<F>& front = L.distFronts[s];
            El::DiagonalSolve
            ( LEFT, NORMAL, front.diag1d, X.distNodes[s-1], true );
        }
    }
}

template<typename F>
inline void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

template<typename F>
inline void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_DIAGONALSOLVE_HPP
