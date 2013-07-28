/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LDL_DISTFRONT_HPP
#define CLIQ_NUMERIC_LDL_DISTFRONT_HPP

namespace cliq {

template<typename F> 
void FrontLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {

template<typename F> 
inline void FrontLDLGeneral
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    CallStackEntry entry("internal::FrontLDLGeneral");
    if( ABR.Height() != ABR.Width() )
        throw std::logic_error("ABR must be square");
    if( AL.Height() != AL.Width()+ABR.Height() )
        throw std::logic_error("AL and ABR must have compatible dimensions");
    if( AL.Grid() != ABR.Grid() )
        throw std::logic_error("AL and ABR must use the same grid");
    if( ABR.ColAlignment() !=
        (AL.ColAlignment()+AL.Width()) % AL.Grid().Height() )
        throw std::logic_error
        ("AL and ABR must have compatible col alignments");
    if( ABR.RowAlignment() != 
        (AL.RowAlignment()+AL.Width()) % AL.Grid().Width() )
        throw std::logic_error
        ("AL and ABR must have compatible row alignments");
    if( orientation == NORMAL )
        throw std::logic_error("FrontLDL must be (conjugate-)transposed.");
#endif
    const Grid& g = AL.Grid();

    // Matrix views
    DistMatrix<F>
        ALTL(g), ALTR(g),  AL00(g), AL01(g), AL02(g),
        ALBL(g), ALBR(g),  AL10(g), AL11(g), AL12(g),
                           AL20(g), AL21(g), AL22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> AL21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21AdjOrTrans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), 
                  AL22B(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
               /**/         AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        AL21_VC_STAR.AlignWith( AL22 );
        AL21_VR_STAR.AlignWith( AL22 );
        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21AdjOrTrans_STAR_MR.AlignWith( AL22 );
        //--------------------------------------------------------------------//
        AL11_STAR_STAR = AL11; 
        elem::LocalLDL( orientation, AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR = AL21;
        elem::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( AL21_VC_STAR );
        elem::DiagonalSolve
        ( RIGHT, NORMAL, d1_STAR_STAR, AL21_VC_STAR );
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
        elem::LocalTrrk
        ( LOWER, orientation, F(-1), leftL, rightL, F(1), AL22T );
        elem::LocalGemm
        ( orientation, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        elem::LocalTrrk( LOWER, orientation, F(-1), leftR, rightR, F(1), ABR );

        elem::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
               /**/         AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
}

template<typename F>
inline void FrontLDLSquare
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    CallStackEntry entry("internal::FrontLDLSquare");
    if( ABR.Height() != ABR.Width() )
        throw std::logic_error("ABR must be square");
    if( AL.Height() != AL.Width()+ABR.Height() )
        throw std::logic_error("AL & ABR must have compatible dimensions");
    if( AL.Grid() != ABR.Grid() )
        throw std::logic_error("AL & ABR must use the same grid");
    if( ABR.ColAlignment() !=
        (AL.ColAlignment()+AL.Width()) % AL.Grid().Height() )
        throw std::logic_error("AL & ABR must have compatible col alignments");
    if( ABR.RowAlignment() != 
        (AL.RowAlignment()+AL.Width()) % AL.Grid().Width() )
        throw std::logic_error("AL & ABR must have compatible row alignments");
    if( orientation == NORMAL )
        throw std::logic_error("FrontLDL must be (conjugate-)transposed.");
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
    DistMatrix<F>
        ALTL(g), ALTR(g),  AL00(g), AL01(g), AL02(g),
        ALBL(g), ALBR(g),  AL10(g), AL11(g), AL12(g),
                           AL20(g), AL21(g), AL22(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> AL11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> AL21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > AL21AdjOrTrans_STAR_MR(g);

    DistMatrix<F,STAR,MC> leftL(g), leftR(g);
    DistMatrix<F,STAR,MR> rightL(g), rightR(g);
    DistMatrix<F> AL22T(g), 
                  AL22B(g);
    DistMatrix<F> AR2T(g), 
                  AR2B(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
               /**/         AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        AL21_VC_STAR.AlignWith( AL22 );
        S21Trans_STAR_MC.AlignWith( AL22 );
        AL21AdjOrTrans_STAR_MR.AlignWith( AL22 );
        //--------------------------------------------------------------------//
        AL11_STAR_STAR = AL11; 
        elem::LocalLDL( orientation, AL11_STAR_STAR, d1_STAR_STAR );
        AL11 = AL11_STAR_STAR;

        AL21_VC_STAR = AL21;
        elem::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          F(1), AL11_STAR_STAR, AL21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( AL21_VC_STAR );
        // SendRecv to form AL21^T[* ,MR] from S21^T[* ,MC], then conjugate
        // if necessary.
        AL21AdjOrTrans_STAR_MR.ResizeTo( AL21.Width(), AL21.Height() );
        {
            if( onDiagonal )
            {
                const int size = AL21.LocalHeight()*AL11.Width();    
                elem::MemCopy
                ( AL21AdjOrTrans_STAR_MR.Buffer(), 
                  S21Trans_STAR_MC.Buffer(), size );
            }
            else
            {
                const int sendSize = AL21.LocalHeight()*AL11.Width();
                const int recvSize = 
                    (AL22.LocalWidth()+ABR.LocalWidth())*AL11.Height();
                // We know that the ldim is the height since we have manually
                // created both temporary matrices.
                mpi::SendRecv
                ( S21Trans_STAR_MC.Buffer(), sendSize, transposeRank, 0,
                  AL21AdjOrTrans_STAR_MR.Buffer(), 
                  recvSize, transposeRank, 0, g.VCComm() );
            }
            elem::DiagonalSolve
            ( LEFT, NORMAL, d1_STAR_STAR, AL21AdjOrTrans_STAR_MR );
            if( orientation == ADJOINT )
                elem::Conjugate( AL21AdjOrTrans_STAR_MR );
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
        elem::LocalTrrk
        ( LOWER, orientation, F(-1), leftL, rightL, F(1), AL22T );
        elem::LocalGemm
        ( orientation, NORMAL, F(-1), leftR, rightL, F(1), AL22B );
        elem::LocalTrrk( LOWER, orientation, F(-1), leftR, rightR, F(1), ABR );

        elem::DiagonalSolve
        ( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        AL21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
               /**/         AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
}

} // namespace internal

template<typename F> 
inline void FrontLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    CallStackEntry entry("FrontLDL");
#endif
    const Grid& grid = AL.Grid();
    if( grid.Height() == grid.Width() )
        internal::FrontLDLSquare( orientation, AL, ABR );
    else
        internal::FrontLDLGeneral( orientation, AL, ABR );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_DISTFRONT_HPP
