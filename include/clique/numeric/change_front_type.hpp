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
    CallStackEntry entry("ChangeFrontType");
#endif
    // Check if this call can be a no-op
    if( frontType == L.frontType ) 
        return;
    const int numDistNodes = L.distFronts.size();    
    DistSymmFront<F>& leafFront = L.distFronts[0];
    const SymmFrontType oldFrontType = L.frontType;

    if( (frontType == LDL_1D        && oldFrontType == LDL_2D)        ||
        (frontType == LDL_SELINV_1D && oldFrontType == LDL_SELINV_2D) ||
        (frontType == SYMM_1D       && oldFrontType == SYMM_2D) )
    {
        // 2d -> 1d
        leafFront.front1dL.LockedAttach
        ( leafFront.front2dL.Height(), 
          leafFront.front2dL.Width(), 0,
          leafFront.front2dL.LockedBuffer(), 
          leafFront.front2dL.LDim(),
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
    else if( (frontType == LDL_2D        && oldFrontType == LDL_1D)        || 
             (frontType == LDL_SELINV_2D && oldFrontType == LDL_SELINV_1D) ||
             (frontType == SYMM_2D       && oldFrontType == SYMM_1D)       || 
             // The following reset the fronts for manual modification
             (frontType == SYMM_2D       && oldFrontType == LDL_1D)        || 
             (frontType == SYMM_2D       && oldFrontType == LDL_SELINV_1D) ||
             (frontType == SYMM_2D       && oldFrontType == BLOCK_LDL_2D) )
    {
        // 1d -> 2d
        leafFront.front2dL.LockedAttach
        ( leafFront.front1dL.Height(), 
          leafFront.front1dL.Width(), 0, 0,
          leafFront.front1dL.LockedBuffer(), 
          leafFront.front1dL.LDim(),
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
            // Invert the unit-diagonal lower triangle
            DistSymmFront<F>& front = L.distFronts[s];
            const int snSize = front.front2dL.Width();
            DistMatrix<F> LT( front.front2dL.Grid() );
            View( LT, front.front2dL, 0, 0, snSize, snSize );
            elem::TriangularInverse( LOWER, UNIT, LT );
        }
    }
    else if( frontType == LDL_SELINV_1D && oldFrontType == LDL_2D )
    {
        // Perform selective inversion and then redistribute to 1d
        leafFront.front1dL.LockedAttach
        ( leafFront.front2dL.Height(), 
          leafFront.front2dL.Width(), 0,
          leafFront.front2dL.LockedBuffer(), 
          leafFront.front2dL.LDim(),
          leafFront.front2dL.Grid() );
        for( int s=1; s<numDistNodes; ++s )
        {
            // Invert the lower triangle
            DistSymmFront<F>& front = L.distFronts[s];
            const Grid& grid = front.front2dL.Grid();
            const int snSize = front.front2dL.Width();
            DistMatrix<F> LT( grid );
            View( LT, front.front2dL, 0, 0, snSize, snSize );
            elem::TriangularInverse( LOWER, UNIT, LT );

            // Copy the data over
            front.front1dL.Empty();
            front.front1dL.SetGrid( grid );
            front.front1dL = front.front2dL;
            front.front2dL.Empty();
        }
    }
    else
        throw std::logic_error("Unavailable front type change");
    L.frontType = frontType;
}

} // namespace cliq
