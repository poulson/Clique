/*
   Modification of include/elemental/advanced/LDL.hpp from Elemental.
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

template<typename F> // F represents a real or complex field
void clique::numeric::LocalFrontLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& AR )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLDL");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("AL and AR must be the same height");
    if( AL.Height() != AL.Width() + AR.Width() )
        throw std::logic_error("[AL AR] must be square");
    if( orientation == NORMAL )
        throw std::logic_error("LocalFrontLDL must be (conjugate-)transposed.");
#endif
    Matrix<F>
        ALTL, ALTR,  AL00, AL01, AL02,
        ALBL, ALBR,  AL10, AL11, AL12,
                     AL20, AL21, AL22;
    Matrix<F>
        ART,  AR0,
        ARB,  AR1,
              AR2;
    Matrix<F> d1;
    Matrix<F> S21;

    Matrix<F> S21T,
              S21B;
    Matrix<F> AL21T,
              AL21B;
    Matrix<F> AR2T,
              AR2B;

    // Start the algorithm
    PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    PartitionDown
    ( AR, ART,
          ARB, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
                /**/        AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        RepartitionDown
        ( ART,  AR0,
         /***/ /***/
                AR1,
          ARB,  AR2, AL11.Height() );

        //--------------------------------------------------------------------//
        // This routine is unblocked, hence the need for us to generalize to 
        // an (ideally) faster blocked algorithm.
        advanced::internal::LDLVar3( orientation, AL11, d1 );

        basic::Trsm( RIGHT, LOWER, orientation, UNIT, (F)1, AL11, AL21 );

        S21 = AL21;
        basic::DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        // For now, perform about 2x as much work as necessary on the 
        // symmetric updates. Eventually, these should be replaced with 
        // custom routines.
        PartitionDown
        ( S21, S21T,
               S21B, AL22.Width() );
        PartitionDown
        ( AL21, AL21T,
                AL21B, AL22.Width() );
        PartitionDown
        ( AR2, AR2T,
               AR2B, AL22.Width() );
        basic::Gemm( NORMAL, orientation, (F)-1, S21, AL21T, (F)1, AL22 );
        basic::Gemm( NORMAL, orientation, (F)-1, S21B, AL21B, (F)1, AR2B );
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( ART,  AR0,
                AR1,
         /***/ /***/
          ARB,  AR2 );

        SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
                /**/        AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<float>& AL, Matrix<float>& AR);

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<double>& AL, Matrix<double>& AR );

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<std::complex<float> >& AL, Matrix<std::complex<float> >& AR );

template void clique::numeric::LocalFrontLDL
( Orientation orientation, 
  Matrix<std::complex<double> >& AL, Matrix<std::complex<double> >& AR );
