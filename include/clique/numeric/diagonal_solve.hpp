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
    const int width = X.Width();
    Matrix<F> XSub;
    for( int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        View( XSub, X.multiVec, node.myOffset, 0, node.size, width );

        Matrix<F> frontTL;
        LockedView( frontTL, frontL, 0, 0, node.size, node.size );
        Matrix<F> d;
        frontTL.GetDiagonal( d );
        elem::DiagonalSolve( LEFT, NORMAL, d, XSub, true );
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
    const int width = X.Width();

    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& front = L.distFronts[s];

        Matrix<F> localXT;
        View
        ( localXT, X.multiVec, 
          node.solveMeta1d.localOffset, 0, node.solveMeta1d.localSize, width );

        elem::DiagonalSolve
        ( LEFT, NORMAL, front.diag.LockedMatrix(), localXT, true );
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
