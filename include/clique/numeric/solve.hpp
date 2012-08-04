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

} // namespace cliq
