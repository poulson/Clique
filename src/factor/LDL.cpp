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

template<typename F>
void clique::factor::LDLT( DistMatrix<F,MC,MR>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("factor::LDLT");
#endif
    clique::factor::LDL( TRANSPOSE, A, supernodeSize );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::factor::LDLH( DistMatrix<F,MC,MR>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("factor::LDLH");
#endif
    clique::factor::LDL( ADJOINT, A, supernodeSize );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::factor::LDL
( Orientation orientation, DistMatrix<F,MC,MR>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("factor::LDL");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed.");
    if( supernodeSize > A.Height() )
        throw std::logic_error("Supernode is too big");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21AdjOrTrans_STAR_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < supernodeSize )
    {
        const int blocksize = std::min(Blocksize(),supernodeSize-A00.Height());
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22, blocksize );

        A21_VC_STAR.AlignWith( A22 );
        A21_VR_STAR.AlignWith( A22 );
        S21Trans_STAR_MC.AlignWith( A22 );
        A21AdjOrTrans_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11; 
        advanced::internal::LocalLDL
        ( orientation, A11_STAR_STAR, d1_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR = A21;
        basic::internal::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          (F)1, A11_STAR_STAR, A21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        basic::DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, A21_VC_STAR );
        A21_VR_STAR = A21_VC_STAR;
        if( orientation == ADJOINT )
            A21AdjOrTrans_STAR_MR.AdjointFrom( A21_VR_STAR );
        else
            A21AdjOrTrans_STAR_MR.TransposeFrom( A21_VR_STAR );
        A12 = A21AdjOrTrans_STAR_MR; // save the (conjugate-)transpose of A21

        basic::internal::LocalTriangularRankK
        ( LOWER, TRANSPOSE,
          (F)-1, S21Trans_STAR_MC, A21AdjOrTrans_STAR_MR, (F)1, A22 );

        basic::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        A21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        A21_VC_STAR.FreeAlignments();
        A21_VR_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        A21AdjOrTrans_STAR_MR.FreeAlignments();

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
