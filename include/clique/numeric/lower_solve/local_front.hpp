/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERSOLVE_LOCALFRONT_HPP
#define CLIQ_NUMERIC_LOWERSOLVE_LOCALFRONT_HPP

namespace cliq {

template<typename F>
void FrontLowerForwardSolve( const Matrix<F>& L, Matrix<F>& X );

template<typename F>
void FrontLowerBackwardSolve
( Orientation orientation, const Matrix<F>& L, Matrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void FrontLowerForwardSolve( const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerForwardSolve");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
#endif
    Matrix<F> LT,
              LB;
    LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    Matrix<F> XT, 
              XB;
    PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), LT, XT, true );
    elem::Gemm( NORMAL, NORMAL, F(-1), LB, XT, F(1), XB );
}

template<typename F>
inline void FrontLowerBackwardSolve
( Orientation orientation, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerBackwardSolve");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( orientation == NORMAL )
        LogicError("This solve must be (conjugate-)transposed");
#endif
    Matrix<F> LT,
              LB;
    LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    Matrix<F> XT,
              XB;
    PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Gemm( orientation, NORMAL, F(-1), LB, XB, F(1), XT );
    elem::Trsm( LEFT, LOWER, orientation, NON_UNIT, F(1), LT, XT, true );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERSOLVE_LOCALFRONT_HPP
