/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
using namespace elemental;

namespace internal {
template<typename F> // represents a real or complex ring
void
SetDiagonalToOne( Matrix<F>& D, int diagOffset )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int diagLength = D.DiagonalLength( diagOffset );
    if( diagOffset >= 0 )
        for( int j=0; j<diagLength; ++j )
            D.Set( j, j+diagOffset, (F)1 );
    else
        for( int j=0; j<diagLength; ++j )
            D.Set( j-diagOffset, j, (F)1 );
#ifndef RELEASE
    PopCallStack();
#endif
}
}

template<typename F>
void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize, 
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerMultiplyNormal");
    if( L.Height() != L.Width() || L.Height() != X.Height() || 
        L.Height() < supernodeSize )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  supernodeSize ~ " << supernodeSize << "\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    Matrix<F>* LMod = const_cast<Matrix<F>*>(&L);
    Matrix<F> LT, LB;
    PartitionDown
    ( *LMod, LT,
             LB, supernodeSize );

    Matrix<F> XT, XB;
    PartitionDown
    ( X, XT,
         XB, supernodeSize );

    basic::Gemm( NORMAL, NORMAL, (F)1, LB, XT, (F)1, XB );

    if( diagOffset == 0 )
    {
        basic::Trmm( LEFT, LOWER, NORMAL, diag, (F)1, LT, XT );
    }
    else
    {
        Matrix<F> d;
        if( diag == UNIT )
        {
            LT.GetDiagonal( d, diagOffset );
            internal::SetDiagonalToOne( LT, diagOffset );
        }
        basic::Trmm( LEFT, LOWER, NORMAL, diag, (F)1, LT, XT );
        if( diag == UNIT )
            LT.SetDiagonal( d, diagOffset );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

// No const since we might modify the diagonal
template<typename F>
void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  int supernodeSize, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerMultiplyTranspose");
    if( L.Height() != L.Width() || L.Height() != X.Height() || 
        L.Height() < supernodeSize )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  supernodeSize ~ " << supernodeSize << "\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transposed");
#endif
    Matrix<F>* LMod = const_cast<Matrix<F>*>(&L);
    Matrix<F> LT, LB;
    PartitionDown
    ( *LMod, LT,
             LB, supernodeSize );

    Matrix<F> XT, XB;
    PartitionDown
    ( X, XT,
         XB, supernodeSize );

    if( diagOffset == 0 )
    {
        basic::Trmm( LEFT, LOWER, orientation, diag, (F)1, LT, XT );
    }
    else
    {
        Matrix<F> d;
        if( diag == UNIT )
        {
            LT.GetDiagonal( d, diagOffset );
            internal::SetDiagonalToOne( LT, diagOffset );
        }
        basic::Trmm( LEFT, LOWER, orientation, diag, (F)1, LT, XT );
        if( diag == UNIT )
            LT.SetDiagonal( d, diagOffset );
    }

    basic::Gemm( orientation, NORMAL, (F)1, LB, XB, (F)1, XT );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<float>& L, Matrix<float>& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<float>& L, Matrix<float>& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<double>& L, Matrix<double>& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<double>& L, Matrix<double>& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
