/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERSOLVE_LOCALFRONTBLOCK_HPP
#define CLIQ_NUMERIC_LOWERSOLVE_LOCALFRONTBLOCK_HPP

namespace cliq {

template<typename F>
void FrontBlockLowerForwardSolve( const Matrix<F>& L, Matrix<F>& X );

template<typename F>
void FrontBlockLowerBackwardSolve
( Orientation orientation, const Matrix<F>& L, Matrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void FrontBlockLowerForwardSolve( const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontBlockLowerForwardSolve");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
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

    // XT := inv(ATL) XT
    Matrix<F> YT( XT );
    elem::Gemm( NORMAL, NORMAL, F(1), LT, YT, F(0), XT );

    // XB := XB - LB XT
    elem::Gemm( NORMAL, NORMAL, F(-1), LB, XT, F(1), XB );
}

template<typename F>
inline void FrontBlockLowerBackwardSolve
( Orientation orientation, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontBlockLowerBackwardSolve");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("This solve must be (conjugate-)transposed");
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

    // YT := LB^[T/H] XB
    Matrix<F> YT;
    elem::Gemm( orientation, NORMAL, F(1), LB, XB, YT );

    // XT := XT - inv(ATL) YT
    elem::Gemm( NORMAL, NORMAL, F(-1), LT, YT, F(1), XT );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERSOLVE_LOCALFRONTBLOCK_HPP
