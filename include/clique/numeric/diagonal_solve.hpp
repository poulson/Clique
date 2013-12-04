/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
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
#ifndef RELEASE
    CallStackEntry cse("LocalDiagonalSolve");
#endif
    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
        elem::DiagonalSolve
        ( LEFT, NORMAL, L.localFronts[s].diag, X.localNodes[s], true );
}

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("LocalDiagonalSolve");
#endif
    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
        elem::DiagonalSolve
        ( LEFT, NORMAL, L.localFronts[s].diag, X.localNodes[s], true );
}

template<typename F> 
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistDiagonalSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmFront<F>& front = L.distFronts[s];
        elem::DiagonalSolve
        ( LEFT, NORMAL, front.diag1d, X.distNodes[s-1], true );
    }
}

template<typename F> 
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistDiagonalSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmFront<F>& front = L.distFronts[s];
        elem::DiagonalSolve
        ( LEFT, NORMAL, front.diag1d, X.distNodes[s-1], true );
    }
}

template<typename F>
inline void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DiagonalSolve");
#endif
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

template<typename F>
inline void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DiagonalSolve");
#endif
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_DIAGONALSOLVE_HPP
