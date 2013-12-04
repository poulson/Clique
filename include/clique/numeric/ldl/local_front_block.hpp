/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
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
void FrontBlockLDL( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void FrontBlockLDL( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry cse("FrontBlockLDL");
#endif
    Matrix<F> ATL,
              ABL;
    PartitionDown
    ( AL, ATL,
          ABL, AL.Width() );
    
    // Make a copy of the original contents of ABL
    Matrix<F> BBL( ABL );

    // Call the standard routine
    FrontLDL( AL, ABR, conjugate );

    // Copy the original contents of ABL back
    ABL = BBL;

    // Overwrite ATL with inv(L D L^[T/H]) = L^[-T/H] D^{-1} L^{-1}
    elem::TriangularInverse( LOWER, UNIT, ATL );
    elem::Trdtrmm( LOWER, ATL, conjugate );
    elem::MakeTrapezoidal( LOWER, ATL );
    Matrix<F> ATLTrans;
    elem::Transpose( ATL, ATLTrans, conjugate );
    elem::MakeTrapezoidal( UPPER, ATLTrans, 1 );
    elem::Axpy( F(1), ATLTrans, ATL );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_LOCALFRONTBLOCK_HPP
