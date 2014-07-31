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
#ifndef CLIQ_NUMERIC_LDL_LOCALFRONTBLOCK_HPP
#define CLIQ_NUMERIC_LDL_LOCALFRONTBLOCK_HPP

namespace cliq {

template<typename F> 
void FrontBlockLDL
( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate=false, bool intraPiv=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void FrontBlockLDL
( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate, bool intraPiv )
{
    DEBUG_ONLY(CallStackEntry cse("FrontBlockLDL"))
    Matrix<F> ATL, ABL;
    PartitionDown( AL, ATL, ABL, AL.Width() );
    
    // Make a copy of the original contents of ABL
    Matrix<F> BBL( ABL );

    if( intraPiv )
    {
        Matrix<Int> p;
        Matrix<F> dSub;
        // TODO: Expose the pivot type as an option?
        El::LDL( ATL, dSub, p, conjugate, El::BUNCH_KAUFMAN_A );

        // Solve against ABL and update ABR
        // NOTE: This does not exploit symmetry
        El::ldl::SolveAfter( ATL, dSub, p, ABL, conjugate );
        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        El::Gemm( NORMAL, orientation, F(-1), ABL, BBL, F(1), ABR );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        El::TriangularInverse( LOWER, UNIT, ATL );
        El::Trdtrmm( LOWER, ATL, dSub, conjugate );
        // TODO: SymmetricPermutation
        El::MakeSymmetric( LOWER, ATL, conjugate );
        El::PermuteRows( ATL, p );
        El::PermuteCols( ATL, p );
    }
    else
    {
        // Call the standard routine
        FrontLDL( AL, ABR, conjugate );

        // Copy the original contents of ABL back
        ABL = BBL;

        // Finish inverting ATL
        El::TriangularInverse( LOWER, UNIT, ATL );
        El::Trdtrmm( LOWER, ATL, conjugate );
        El::MakeSymmetric( LOWER, ATL, conjugate );
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_LOCALFRONTBLOCK_HPP
