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
void ChangeFrontType( DistSymmFrontTree<F>& L, SymmFrontType frontType );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

// This routine could be modified later so that it uses much less memory
// by replacing the '=' redistributions with piece-by-piece redistributions.
template<typename F>
inline void ChangeFrontType( DistSymmFrontTree<F>& L, SymmFrontType frontType )
{
#ifndef RELEASE
    PushCallStack("ChangeFrontType");
#endif
    // Check if this call can be a no-op
    if( frontType == L.frontType ) 
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }
    const int numDistNodes = L.distFronts.size();    
    DistSymmFront<F>& leafFront = L.distFronts[0];
    const SymmFrontType oldFrontType = L.frontType;

    if( (frontType == LDL_1D && oldFrontType == LDL_2D) ||
        (frontType == LDL_SELINV_1D && oldFrontType == LDL_SELINV_2D) ||
        (frontType == SYMM_1D && oldFrontType == SYMM_2D) )

    {
        // 2d -> 1d
        leafFront.front1dL.LockedView
        ( leafFront.front2dL.Height(), 
          leafFront.front2dL.Width(), 0,
          leafFront.front2dL.LockedLocalBuffer(), 
          leafFront.front2dL.LocalLDim(),
          leafFront.front2dL.Grid() );
        for( int s=1; s<numDistNodes; ++s )
        {
            DistSymmFront<F>& front = L.distFronts[s];
            front.front1dL.Empty();
            front.front1dL.SetGrid( front.front2dL.Grid() );
            front.front1dL = front.front2dL;
            front.front2dL.Empty();
        }
    }
    else if( (frontType == LDL_2D && oldFrontType == LDL_1D) || 
             (frontType == LDL_SELINV_2D && oldFrontType == LDL_SELINV_1D) ||
             (frontType == SYMM_2D && oldFrontType == SYMM_1D) )
    {
        // 1d -> 2d
        leafFront.front2dL.LockedView
        ( leafFront.front1dL.Height(), 
          leafFront.front1dL.Width(), 0, 0,
          leafFront.front1dL.LockedLocalBuffer(), 
          leafFront.front1dL.LocalLDim(),
          leafFront.front1dL.Grid() );
        for( int s=1; s<numDistNodes; ++s )
        {
            DistSymmFront<F>& front = L.distFronts[s];
            front.front2dL.Empty();
            front.front2dL.SetGrid( front.front1dL.Grid() );
            front.front2dL = front.front1dL;
            front.front1dL.Empty();
        }
    }
    else if( frontType == LDL_SELINV_2D && oldFrontType == LDL_2D )
    {
        // Perform selective inversion
        for( int s=1; s<numDistNodes; ++s )
        {
            DistSymmFront<F>& front = L.distFronts[s];
            const Grid& grid = front.front2dL.Grid();
            const int snSize = front.front2dL.Width();

            // Invert the strictly lower portion of the diagonal block, and
            // then make the strictly upper triangle zero
            DistMatrix<F> LT( grid );
            LT.View( front.front2dL, 0, 0, snSize, snSize );
            elem::TriangularInverse( LOWER, UNIT, LT );
            elem::MakeTrapezoidal( LEFT, LOWER, 0, LT );
        }
    }
    else if( frontType == LDL_SELINV_1D && oldFrontType == LDL_2D )
    {
        // Perform selective inversion and then redistribute to 1d
        leafFront.front1dL.LockedView
        ( leafFront.front2dL.Height(), 
          leafFront.front2dL.Width(), 0,
          leafFront.front2dL.LockedLocalBuffer(), 
          leafFront.front2dL.LocalLDim(),
          leafFront.front2dL.Grid() );
        for( int s=1; s<numDistNodes; ++s )
        {
            DistSymmFront<F>& front = L.distFronts[s];
            const Grid& grid = front.front2dL.Grid();
            const int snSize = front.front2dL.Width();

            // Invert the strictly lower portion of the diagonal block 
            DistMatrix<F> LT( grid );
            LT.View( front.front2dL, 0, 0, snSize, snSize );
            elem::TriangularInverse( LOWER, UNIT, LT );

            // Copy the data and make the strictly upper triangle zero
            front.front1dL.Empty();
            front.front1dL.SetGrid( grid );
            front.front1dL = front.front2dL;
            front.front2dL.Empty();
            elem::MakeTrapezoidal( LEFT, LOWER, 0, front.front1dL );
        }
    }
    else
        throw std::logic_error("Unavailable front type change");
    L.frontType = frontType;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
