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
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X );

template<typename F>
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X );

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("LocalDiagonalSolve");
#endif
    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;

        Matrix<F> frontTL;
        LockedView( frontTL, frontL, 0, 0, node.size, node.size );
        Matrix<F> d;
        frontTL.GetDiagonal( d );
        elem::DiagonalSolve( LEFT, NORMAL, d, X.localNodes[s], true );
    }
}

template<typename F> 
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistDiagonalSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmFront<F>& front = L.distFronts[s];
        elem::DiagonalSolve
        ( LEFT, NORMAL, front.diag.LockedMatrix(), X.distNodes[s-1], true );
    }
}

template<typename F>
inline void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, 
  DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DiagonalSolve");
#endif
    LocalDiagonalSolve( info, L, X );
    DistDiagonalSolve( info, L, X );
}

} // namespace cliq
