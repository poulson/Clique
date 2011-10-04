/*
   Modification of include/elemental/basic/level3/Trsm/TrsmLLN.hpp 
   from Elemental.
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

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

template<typename F>
void clique::numeric::LocalFrontLDLDiagonalSolve
( int supernodeSize,
  const Matrix<F>& d, Matrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLDLDiagonalSolve");
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    if( d.Height() != X.Height() )
        throw std::logic_error("Invalid height of X");
#endif
    basic::DiagonalSolve( LEFT, NORMAL, d, X, checkIfSingular );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::LocalFrontLDLForwardSolve
( int supernodeSize, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLDLForwardSolve");
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
#endif
    // Matrix views
    Matrix<F>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

    // Start the algorithm
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XT.Height() < supernodeSize )
    {
        const int blocksize = std::min(Blocksize(),supernodeSize-XT.Height());
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22, blocksize );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, blocksize );

        //--------------------------------------------------------------------//
        // X1 := L11^-1 X1
        basic::Trsm( LEFT, LOWER, NORMAL, UNIT, (F)1, L11, X1 );

        // X2 -= L21 X1
        basic::Gemm( NORMAL, NORMAL, (F)-1, L21, X1, (F)1, X2 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

// The unit upper triangle and its (conjugate-)transpose 
// (with the exception of the diagonal) must be explicitly stored
template<typename F>
void clique::numeric::LocalFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize, const Matrix<F>& U, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLDLBackwardSolve");
    if( U.Height() != U.Width() || U.Height() != X.Height() || 
        U.Height() < supernodeSize )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  supernodeSize ~ " << supernodeSize << "\n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    // Matrix views
    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, U.Height()-supernodeSize );
    PartitionUp
    ( X, XT,
         XB, X.Height()-supernodeSize );

    // Subtract off the parent updates
    basic::Gemm( NORMAL, NORMAL, (F)-1, UTR, XB, (F)1, XT );

    // Solve the remaining triangular system
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        // X1 := U11^-1 X1
        basic::Trsm( LEFT, UPPER, NORMAL, UNIT, (F)1, U11, X1 );

        // X0 -= U01 X1
        basic::Gemm( NORMAL, NORMAL, (F)-1, U01, X1, (F)1, X0 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLDLForwardSolve
( int supernodeSize,
  const Matrix<float>& L, Matrix<float>& X );
template void clique::numeric::LocalFrontLDLDiagonalSolve
( int supernodeSize,
  const Matrix<float>& d, Matrix<float>& X,
  bool checkIfSingular );
template void clique::numeric::LocalFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const Matrix<float>& U, Matrix<float>& X );

template void clique::numeric::LocalFrontLDLForwardSolve
( int supernodeSize,
  const Matrix<double>& L, Matrix<double>& X );
template void clique::numeric::LocalFrontLDLDiagonalSolve
( int supernodeSize,
  const Matrix<double>& d, Matrix<double>& X,
  bool checkIfSingular );
template void clique::numeric::LocalFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const Matrix<double>& U, Matrix<double>& X );

template void clique::numeric::LocalFrontLDLForwardSolve
( int supernodeSize,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalFrontLDLDiagonalSolve
( int supernodeSize,
  const Matrix<std::complex<float> >& d, Matrix<std::complex<float> >& X,
  bool checkIfSingular );
template void clique::numeric::LocalFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const Matrix<std::complex<float> >& U, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalFrontLDLForwardSolve
( int supernodeSize,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalFrontLDLDiagonalSolve
( int supernodeSize,
  const Matrix<std::complex<double> >& d, Matrix<std::complex<double> >& X,
  bool checkIfSingular );
template void clique::numeric::LocalFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const Matrix<std::complex<double> >& U, Matrix<std::complex<double> >& X );
