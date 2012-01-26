/*
   Modification of include/elemental/basic/level3/Trsm/TrsmLLN.hpp 
   from Elemental.
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

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
#include "clique.hpp"

namespace cliq {

template<typename F>
void numeric::LocalFrontLowerForwardSolve
( Diagonal diag, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalFrontLowerForwardSolve");
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
    elem::LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    Matrix<F> XT, 
              XB;
    elem::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Trsm( LEFT, LOWER, NORMAL, diag, (F)1, LT, XT, true );
    elem::Gemm( NORMAL, NORMAL, (F)-1, LB, XT, (F)1, XB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag, 
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalFrontLowerBackwardSolve");
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
    elem::LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    Matrix<F> XT,
              XB;
    elem::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Gemm( orientation, NORMAL, (F)-1, LB, XB, (F)1, XT );
    elem::Trsm( LEFT, LOWER, orientation, diag, (F)1, LT, XT, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

template void cliq::numeric::LocalFrontLowerForwardSolve
( Diagonal diag, 
  const Matrix<float>& L, Matrix<float>& X );
template void cliq::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<float>& L, Matrix<float>& X );

template void cliq::numeric::LocalFrontLowerForwardSolve
( Diagonal diag, 
  const Matrix<double>& L, Matrix<double>& X );
template void cliq::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<double>& L, Matrix<double>& X );

template void cliq::numeric::LocalFrontLowerForwardSolve
( Diagonal diag,
  const Matrix<Complex<float> >& L, Matrix<Complex<float> >& X );
template void cliq::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<Complex<float> >& L, Matrix<Complex<float> >& X );

template void cliq::numeric::LocalFrontLowerForwardSolve
( Diagonal diag,
  const Matrix<Complex<double> >& L, Matrix<Complex<double> >& X );
template void cliq::numeric::LocalFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const Matrix<Complex<double> >& L, Matrix<Complex<double> >& X );
