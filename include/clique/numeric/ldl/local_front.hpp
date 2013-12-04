/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LDL_LOCALFRONT_HPP
#define CLIQ_NUMERIC_LDL_LOCALFRONT_HPP

namespace cliq {

template<typename F> 
void FrontLDL( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void FrontLDL( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLDL");
    if( ABR.Height() != ABR.Width() )
        LogicError("ABR must be square");
    if( AL.Height() != AL.Width() + ABR.Width() )
        LogicError("AL and ABR don't have conformal dimensions");
#endif
    Matrix<F>
        ALTL, ALTR,  AL00, AL01, AL02,
        ALBL, ALBR,  AL10, AL11, AL12,
                     AL20, AL21, AL22;
    Matrix<F> d1;
    Matrix<F> S21;

    Matrix<F> S21T,
              S21B;
    Matrix<F> AL21T,
              AL21B;

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Start the algorithm
    PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
                /**/        AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        //--------------------------------------------------------------------//
        // This routine is unblocked, hence the need for us to generalize to 
        // an (ideally) faster blocked algorithm.
        elem::ldl::Var3( AL11, conjugate );
        AL11.GetDiagonal( d1 );

        elem::Trsm( RIGHT, LOWER, orientation, UNIT, F(1), AL11, AL21 );

        S21 = AL21;
        elem::DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        PartitionDown
        ( S21, S21T,
               S21B, AL22.Width() );
        PartitionDown
        ( AL21, AL21T,
                AL21B, AL22.Width() );
        elem::Gemm( NORMAL, orientation, F(-1), S21, AL21T, F(1), AL22 );
        elem::MakeTriangular( LOWER, AL22 );
        elem::internal::TrrkNT
        ( LOWER, orientation, F(-1), S21B, AL21B, F(1), ABR );
        //--------------------------------------------------------------------//

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
                /**/        AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_LOCALFRONT_HPP
