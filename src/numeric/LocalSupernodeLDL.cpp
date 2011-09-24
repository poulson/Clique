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
void clique::numeric::LocalSupernodeLDL
( Orientation orientation, Matrix<F>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalSupernodeLDL");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( orientation == NORMAL )
        throw std::logic_error
        ("LocalSupernodeLDL must be (conjugate-)transposed.");
    if( supernodeSize > A.Height() )
        throw std::logic_error("Supernode is too big");
#endif
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<F> d1;
    Matrix<F> S21;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < supernodeSize )
    {
        const int blocksize = std::min(Blocksize(),supernodeSize-ATL.Height());
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22, blocksize );

        //--------------------------------------------------------------------//
        // This routine is unblocked, hence the need for us to generalize to 
        // an (ideally) faster blocked algorithm.
        advanced::internal::LDLVar3( orientation, A11, d1 );

        basic::Trsm( RIGHT, LOWER, orientation, UNIT, (F)1, A11, A21 );

        S21 = A21;
        basic::DiagonalSolve( RIGHT, NORMAL, d1, A21 );
        basic::Transpose( A21, A12 );

        // For now, just perform 2x as much work as necessary via a gemm.
        // Eventually, this should be replaced with a custom routine.
        basic::Gemm( NORMAL, orientation, (F)-1, S21, A21, (F)1, A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalSupernodeLDL
( Orientation orientation, Matrix<float>& A, int supernodeSize );

template void clique::numeric::LocalSupernodeLDL
( Orientation orientation, Matrix<double>& A, int supernodeSize );

template void clique::numeric::LocalSupernodeLDL
( Orientation orientation, 
  Matrix<std::complex<float> >& A, int supernodeSize );

template void clique::numeric::LocalSupernodeLDL
( Orientation orientation, 
  Matrix<std::complex<double> >& A, int supernodeSize );
