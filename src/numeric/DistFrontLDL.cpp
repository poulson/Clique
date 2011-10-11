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
void clique::numeric::DistFrontLDL
( Orientation orientation, DistMatrix<F,MC,MR>& AL, DistMatrix<F,MC,MR>& AR )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontLDL");
#endif
    const Grid& grid = AL.Grid();
    if( grid.Height() == grid.Width() )
        internal::DistFrontLDLSquare( orientation, AL, AR );
    else
        internal::DistFrontLDLGeneral( orientation, AL, AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::internal::DistFrontLDLGeneral
( Orientation orientation, DistMatrix<F,MC,MR>& AL, DistMatrix<F,MC,MR>& AR )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::internal::DistFrontLDLGeneral");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("AL and AR must be the same height");
    if( AL.Height() != AL.Width() + AR.Width() )
        throw std::logic_error
        ("The combined widths of AL and AR must equal their heights");
    if( AL.Grid() != AR.Grid() )
        throw std::logic_error("AL and AR must use the same grid");
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw std::logic_error("AL and AR must have the same col alignments");
    if( AR.RowAlignment() != 
        (AL.RowAlignment()+AL.Width()) % AL.Grid().Width() )
        throw std::logic_error("AL and AR must have compatible row alignments");
    if( orientation == NORMAL )
        throw std::logic_error("DistFrontLDL must be (conjugate-)transposed.");
#endif
    const Grid& g = AL.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ALTL(g), ALTR(g),  AL00(g), AL01(g), AL02(g),
        ALBL(g), ALBR(g),  AL10(g), AL11(g), AL12(g),
                           AL20(g), AL21(g), AL22(g);
    DistMatrix<F,MC,MR>
        ART(g),  AR0(g),
        ARB(g),  AR1(g),
                 AR2(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> AL21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21AdjOrTrans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F,MC,  MR> AL22T(g), 
                          AL22B(g);
    DistMatrix<F,MC,  MR> AR2T(g), 
                          AR2B(g);

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
               /**/         AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        RepartitionDown
        ( ART,  AR0,
         /***/ /***/
                AR1,
          ARB,  AR2, AL11.Height() );

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VR_STAR.AlignWith( AL22 );
        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21AdjOrTrans_STAR_MR.AlignWith( AL22 );
        //--------------------------------------------------------------------//
        AL11_STAR_STAR = AL11; 
        elemental::advanced::internal::LocalLDL
        ( orientation, AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR = AL21;
        basic::internal::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          (F)1, AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( AL21_VC_STAR );
        basic::DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, AL21_VC_STAR );
        AL21_VR_STAR = AL21_VC_STAR;
        if( orientation == ADJOINT )
            AL21AdjOrTrans_STAR_MR.AdjointFrom( AL21_VR_STAR );
        else
            AL21AdjOrTrans_STAR_MR.TransposeFrom( AL21_VR_STAR );

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight
        ( S21Trans_STAR_MC, 
          leftL, leftR, AL22.Width() );
        PartitionRight
        ( AL21AdjOrTrans_STAR_MR,
          rightL, rightR, AL22.Width() );
        PartitionDown
        ( AL22, AL22T,
                AL22B, AL22.Width() );
        PartitionDown
        ( AR2,  AR2T,
                AR2B, AL22.Width() );
        basic::internal::LocalTriangularRankK
        ( LOWER, orientation, 
          (F)-1, leftL, rightL, (F)1, AL22T );
        basic::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, leftR, rightL, (F)1, AL22B );
        basic::internal::LocalTriangularRankK
        ( LOWER, orientation,
          (F)-1, leftR, rightR, (F)1, AR2B );

        basic::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        AL21_VC_STAR.FreeAlignments();
        AL21_VR_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        AL21AdjOrTrans_STAR_MR.FreeAlignments();

        SlidePartitionDown
        ( ART,   AR0,
                 AR1, 
         /***/  /***/
          ARB,   AR2 );

        SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
               /**/         AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::internal::DistFrontLDLSquare
( Orientation orientation, DistMatrix<F,MC,MR>& AL, DistMatrix<F,MC,MR>& AR )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::internal::DistFrontLDLSquare");
    if( AL.Height() != AR.Height() )
        throw std::logic_error("AL and AR must be the same height");
    if( AL.Height() != AL.Width() + AR.Width() )
        throw std::logic_error
        ("The combined widths of AL and AR must equal their heights");
    if( AL.Grid() != AR.Grid() )
        throw std::logic_error("AL and AR must use the same grid");
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw std::logic_error("AL and AR must have the same col alignments");
    if( AR.RowAlignment() != 
        (AL.RowAlignment()+AL.Width()) % AL.Grid().Width() )
        throw std::logic_error("AL and AR must have compatible row alignments");
    if( orientation == NORMAL )
        throw std::logic_error("DistFrontLDL must be (conjugate-)transposed.");
#endif
    const Grid& g = AL.Grid();
#ifndef RELEASE
    if( g.Height() != g.Width() )
        throw std::logic_error("This routine assumes a square process grid");
#endif

    // Find the process holding our transposed data
    const int r = g.Height();
    int transposeRank;
    {
        const int colAlignment = AL.ColAlignment();
        const int rowAlignment = AL.RowAlignment();
        const int colShift = AL.ColShift();
        const int rowShift = AL.RowShift();

        const int transposeRow = (colAlignment+rowShift) % r;
        const int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    // Matrix views
    DistMatrix<F,MC,MR>
        ALTL(g), ALTR(g),  AL00(g), AL01(g), AL02(g),
        ALBL(g), ALBR(g),  AL10(g), AL11(g), AL12(g),
                           AL20(g), AL21(g), AL22(g);
    DistMatrix<F,MC,MR>
        ART(g),  AR0(g),
        ARB(g),  AR1(g),
                 AR2(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21AdjOrTrans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F,MC,  MR> AL22T(g), 
                          AL22B(g);
    DistMatrix<F,MC,  MR> AR2T(g), 
                          AR2B(g);

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
               /**/         AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        RepartitionDown
        ( ART,  AR0,
         /***/ /***/
                AR1,
          ARB,  AR2, AL11.Height() );

        AL21_VC_STAR.AlignWith( AL22 );
        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21AdjOrTrans_STAR_MR.AlignWith( AL22 );
        //--------------------------------------------------------------------//
        AL11_STAR_STAR = AL11; 
        elemental::advanced::internal::LocalLDL
        ( orientation, AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR = AL21;
        basic::internal::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          (F)1, AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( AL21_VC_STAR );
        // SendRecv to form AL21^T[* ,MR] from S21^T[* ,MC], then conjugate
        // if necessary.
        AL21AdjOrTrans_STAR_MR.ResizeTo( AL21.Width(), AL21.Height() );
        {
            if( onDiagonal )
            {
                const int size = AL21.LocalHeight()*AL11.Width();    
                memcpy
                ( AL21AdjOrTrans_STAR_MR.LocalBuffer(), 
                  S21Trans_STAR_MC.LocalBuffer(), size*sizeof(F) );
            }
            else
            {
                const int sendSize = AL21.LocalHeight()*AL11.Width();
                const int recvSize = 
                    (AL22.LocalWidth()+AR2.LocalWidth())*AL11.Height();
                // We know that the ldim is the height since we have manually
                // created both temporary matrices.
                mpi::SendRecv
                ( S21Trans_STAR_MC.LocalBuffer(), sendSize, transposeRank, 0,
                  AL21AdjOrTrans_STAR_MR.LocalBuffer(), 
                  recvSize, transposeRank, 0, g.VCComm() );
            }
            basic::DiagonalSolve
            ( LEFT, NORMAL, d1_STAR_STAR, AL21AdjOrTrans_STAR_MR );
            if( orientation == ADJOINT )
                basic::Conjugate( AL21AdjOrTrans_STAR_MR );
        }

        // Partition the update of the bottom-right corner into three pieces
        PartitionRight
        ( S21Trans_STAR_MC, 
          leftL, leftR, AL22.Width() );
        PartitionRight
        ( AL21AdjOrTrans_STAR_MR,
          rightL, rightR, AL22.Width() );
        PartitionDown
        ( AL22, AL22T,
                AL22B, AL22.Width() );
        PartitionDown
        ( AR2,  AR2T,
                AR2B, AL22.Width() );
        basic::internal::LocalTriangularRankK
        ( LOWER, orientation, 
          (F)-1, leftL, rightL, (F)1, AL22T );
        basic::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, leftR, rightL, (F)1, AL22B );
        basic::internal::LocalTriangularRankK
        ( LOWER, orientation,
          (F)-1, leftR, rightR, (F)1, AR2B );

        basic::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        AL21_VC_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        AL21AdjOrTrans_STAR_MR.FreeAlignments();

        SlidePartitionDown
        ( ART,   AR0,
                 AR1, 
         /***/  /***/
          ARB,   AR2 );

        SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
               /**/         AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::DistFrontLDL
( Orientation orientation, 
  DistMatrix<float,MC,MR>& AL, 
  DistMatrix<float,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLGeneral
( Orientation orientation, 
  DistMatrix<float,MC,MR>& AL, 
  DistMatrix<float,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLSquare
( Orientation orientation, 
  DistMatrix<float,MC,MR>& AL, 
  DistMatrix<float,MC,MR>& AR );

template void clique::numeric::DistFrontLDL
( Orientation orientation, 
  DistMatrix<double,MC,MR>& AL, 
  DistMatrix<double,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLGeneral
( Orientation orientation, 
  DistMatrix<double,MC,MR>& AL, 
  DistMatrix<double,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLSquare
( Orientation orientation, 
  DistMatrix<double,MC,MR>& AL, 
  DistMatrix<double,MC,MR>& AR );

template void clique::numeric::DistFrontLDL
( Orientation orientation, 
  DistMatrix<std::complex<float>,MC,MR>& AL, 
  DistMatrix<std::complex<float>,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLGeneral
( Orientation orientation, 
  DistMatrix<std::complex<float>,MC,MR>& AL, 
  DistMatrix<std::complex<float>,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLSquare
( Orientation orientation, 
  DistMatrix<std::complex<float>,MC,MR>& AL, 
  DistMatrix<std::complex<float>,MC,MR>& AR );

template void clique::numeric::DistFrontLDL
( Orientation orientation, 
  DistMatrix<std::complex<double>,MC,MR>& AL,
  DistMatrix<std::complex<double>,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLGeneral
( Orientation orientation, 
  DistMatrix<std::complex<double>,MC,MR>& AL,
  DistMatrix<std::complex<double>,MC,MR>& AR );
template void clique::numeric::internal::DistFrontLDLSquare
( Orientation orientation, 
  DistMatrix<std::complex<double>,MC,MR>& AL,
  DistMatrix<std::complex<double>,MC,MR>& AR );
