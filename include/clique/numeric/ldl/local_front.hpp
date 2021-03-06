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
#ifndef CLIQ_NUMERIC_LDL_LOCALFRONT_HPP
#define CLIQ_NUMERIC_LDL_LOCALFRONT_HPP

namespace cliq {

template<typename F> 
void FrontLDL( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate=false );

template<typename F>
void FrontLDLIntraPiv
( Matrix<F>& AL, Matrix<F>& subdiag, Matrix<Int>& piv, Matrix<F>& ABR, 
  bool conjugate=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void FrontLDL( Matrix<F>& AL, Matrix<F>& ABR, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLDL");
        if( ABR.Height() != ABR.Width() )
            LogicError("ABR must be square");
        if( AL.Height() != AL.Width() + ABR.Width() )
            LogicError("AL and ABR don't have conformal dimensions");
    )
    const Int m = AL.Height();
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> d1;
    Matrix<F> S21;

    Matrix<F> S21T, S21B;
    Matrix<F> AL21T, AL21B;

    const Int bsize = El::Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = El::Min(bsize,n-k);
        auto AL11 = ViewRange( AL, k,    k,    k+nb, k+nb );
        auto AL21 = ViewRange( AL, k+nb, k,    m,    k+nb );
        auto AL22 = ViewRange( AL, k+nb, k+nb, m,    n    );

        El::LDL( AL11, conjugate );
        AL11.GetDiagonal( d1 );

        El::Trsm( RIGHT, LOWER, orientation, UNIT, F(1), AL11, AL21 );

        S21 = AL21;
        El::DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        PartitionDown( S21, S21T, S21B, AL22.Width() );
        PartitionDown( AL21, AL21T, AL21B, AL22.Width() );
        El::Gemm( NORMAL, orientation, F(-1), S21, AL21T, F(1), AL22 );
        El::MakeTriangular( LOWER, AL22 );
        El::Trrk( LOWER, NORMAL, orientation, F(-1), S21B, AL21B, F(1), ABR );
    }
}

template<typename F>
void FrontLDLIntraPiv
( Matrix<F>& AL, Matrix<F>& subdiag, Matrix<Int>& piv, Matrix<F>& ABR, 
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("FrontLDLIntraPiv"))
    const Int n = AL.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> ATL, ABL;
    PartitionDown( AL, ATL, ABL, n );

    El::LDL( ATL, subdiag, piv, conjugate, El::BUNCH_KAUFMAN_A );
    auto diag = ATL.GetDiagonal();

    El::PermuteCols( ABL, piv );
    El::Trsm( RIGHT, LOWER, orientation, UNIT, F(1), ATL, ABL );
    Matrix<F> SBL( ABL );

    El::QuasiDiagonalSolve( RIGHT, LOWER, diag, subdiag, ABL, conjugate );
    El::Trrk( LOWER, NORMAL, orientation, F(-1), SBL, ABL, F(1), ABR );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_LOCALFRONT_HPP
