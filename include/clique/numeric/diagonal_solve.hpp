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
#ifndef CLIQUE_DIAGONAL_SOLVE_HPP
#define CLIQUE_DIAGONAL_SOLVE_HPP 1

namespace cliq {

template<typename F>
void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

template<typename F>
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

template<typename F>
inline void LocalDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("LocalDiagonalSolve");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    Matrix<F> XSub;
    for( int s=0; s<numLocalNodes; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        XSub.View( X, node.myOffset, 0, node.size, width );

        Matrix<F> frontTL;
        frontTL.LockedView( frontL, 0, 0, node.size, node.size );
        Matrix<F> d;
        frontTL.GetDiagonal( d );
        elem::DiagonalSolve( LEFT, NORMAL, d, XSub, true );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
void DistDiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("DistDiagonalSolve");
#endif
    const int numDistNodes = info.distNodes.size();
    const int width = localX.Width();

    for( int s=1; s<numDistNodes; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const DistSymmFront<F>& front = L.distFronts[s];

        Matrix<F> localXT;
        localXT.View( localX, node.localOffset1d, 0, node.localSize1d, width );

        elem::DiagonalSolve
        ( LEFT, NORMAL, front.diag.LockedLocalMatrix(), localXT, true );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void DiagonalSolve
( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("DiagonalSolve");
#endif
    LocalDiagonalSolve( info, L, localX );
    DistDiagonalSolve( info, L, localX );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // CLIQUE_DIAGONAL_SOLVE_HPP 
