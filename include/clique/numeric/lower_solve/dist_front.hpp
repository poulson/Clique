/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERSOLVE_DISTFRONT_HPP
#define CLIQ_NUMERIC_LOWERSOLVE_DISTFRONT_HPP

namespace cliq {

template<typename F>
void FrontLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool singleL11AllGather=true );
template<typename F>
void FrontLowerForwardSolve( const DistMatrix<F>& L, DistMatrix<F>& X );

template<typename F>
void FrontLowerBackwardSolve
( Orientation orientation, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool singleL11AllGather=true );
template<typename F>
void FrontLowerBackwardSolve
( Orientation orientation, const DistMatrix<F>& L, DistMatrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename F>
inline void ForwardMany
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        elem::LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, 
          F(1), L11_STAR_STAR, X1_STAR_STAR, true );
        X1 = X1_STAR_STAR;

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        elem::LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
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
}

template<typename F>
void FormDiagonalBlocks
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,STAR,STAR>& D, bool conjugate )
{
    const Grid& g = L.Grid();

    const Int height = L.Width();
    const Int blocksize = Blocksize();

    const int commRank = g.VCRank();
    const int commSize = g.Size();

    const Int localHeight = Length(height,commRank,commSize);
    const Int maxLocalHeight = MaxLength(height,commSize);
    const Int portionSize = maxLocalHeight*blocksize;

    std::vector<F> sendBuffer( portionSize );
    const Int colShift = L.ColShift();
    const Int LLDim = L.LDim();
    const F* LBuffer = L.LockedBuffer();
    if( conjugate )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*commSize;
            const Int block = i / blocksize;
            const Int jStart = block*blocksize;
            const Int b = std::min(height-jStart,blocksize);
            for( Int jOff=0; jOff<b; ++jOff )
                sendBuffer[iLoc*blocksize+jOff] = 
                    Conj(LBuffer[iLoc+(jStart+jOff)*LLDim]);
        }
    }
    else
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*commSize;
            const Int block = i / blocksize;
            const Int jStart = block*blocksize;
            const Int b = std::min(height-jStart,blocksize);
            for( Int jOff=0; jOff<b; ++jOff )
                sendBuffer[iLoc*blocksize+jOff] = 
                    LBuffer[iLoc+(jStart+jOff)*LLDim];
        }
    }

    std::vector<F> recvBuffer( portionSize*commSize );
    mpi::AllGather
    ( &sendBuffer[0], portionSize, &recvBuffer[0], portionSize, g.VCComm() );
    std::vector<F>().swap( sendBuffer );
    
    D.ResizeTo( blocksize, height );
    F* DBuffer = D.Buffer();
    const Int DLDim = D.LDim();
    for( Int proc=0; proc<commSize; ++proc )
    {
        const F* procRecv = &recvBuffer[proc*portionSize];
        const Int procLocalHeight = Length(height,proc,commSize);
        for( Int iLoc=0; iLoc<procLocalHeight; ++iLoc )
        {
            const Int i = proc + iLoc*commSize;
            for( Int jOff=0; jOff<blocksize; ++jOff )
                DBuffer[jOff+i*DLDim] = procRecv[jOff+iLoc*blocksize];
        }
    }
}

template<typename F>
void AccumulateRHS( const DistMatrix<F,VC,STAR>& X, DistMatrix<F,STAR,STAR>& Z )
{
    const Int height = X.Height();
    const Int width = X.Width();
    Z.Empty();
    elem::Zeros( Z, height, width );

    const Int localHeight = X.LocalHeight();
    const Int colShift = X.ColShift();
    const int commSize = X.Grid().Size();
    const F* XBuffer = X.LockedBuffer();
    F* ZBuffer = Z.Buffer();
    const Int XLDim = X.LDim();
    const Int ZLDim = Z.LDim();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = colShift + iLoc*commSize;
        for( Int j=0; j<width; ++j )
            ZBuffer[i+j*ZLDim] = XBuffer[iLoc+j*XLDim];
    }
    mpi::AllReduce( ZBuffer, ZLDim*width, mpi::SUM, X.Grid().VCComm() );
}

template<typename F>
void ForwardSingle( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11Trans_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    DistMatrix<F,STAR,STAR> D(g);
    FormDiagonalBlocks( L, D, false );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        LockedView
        ( L11Trans_STAR_STAR, D, 0, L00.Height(), L11.Height(), L11.Height() );

        //--------------------------------------------------------------------//
        AccumulateRHS( X1, X1_STAR_STAR ); // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        elem::LocalTrsm
        ( LEFT, UPPER, TRANSPOSE, NON_UNIT, 
          F(1), L11Trans_STAR_STAR, X1_STAR_STAR, true );
        X1 = X1_STAR_STAR;

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        elem::LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
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
}

template<typename F>
void BackwardMany
( Orientation orientation, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    // TODO: Replace this with modified inline code?
    elem::internal::TrsmLLTSmall( orientation, NON_UNIT, F(1), L, X, true );
}

template<typename F>
void BackwardSingle
( Orientation orientation,
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> L11AdjOrTrans_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> Z1_STAR_STAR(g);

    DistMatrix<F,STAR,STAR> D(g);
    if( orientation == TRANSPOSE ) 
        FormDiagonalBlocks( L, D, false );
    else 
        FormDiagonalBlocks( L, D, true );
    const Int blocksize = Blocksize();
    const Int firstBlocksize = 
        ( L.Height()%blocksize==0 ?
          blocksize :
          L.Height()%blocksize );

    // Start the algorithm
    Int b = firstBlocksize;
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22, b );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2, b );

        LockedView( L11AdjOrTrans_STAR_STAR, D, 0, L00.Height(), b, b );

        //--------------------------------------------------------------------//
        // X1 -= L21' X2
        elem::LocalGemm( orientation, NORMAL, F(-1), L21, X2, Z1_STAR_STAR );
        elem::internal::AddInLocalData( X1, Z1_STAR_STAR );
        Z1_STAR_STAR.SumOverGrid();

        // X1 := L11^-1 X1
        elem::LocalTrsm
        ( LEFT, UPPER, NORMAL, UNIT, 
          F(1), L11AdjOrTrans_STAR_STAR, Z1_STAR_STAR );
        X1 = Z1_STAR_STAR;
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        b = blocksize;
    }
}

} // namespace internal

template<typename F>
inline void FrontLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool singleL11AllGather )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerForwardSolve");
    if( L.Grid() != X.Grid() )
        LogicError("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( L.ColAlign() != X.ColAlign() )
        LogicError("L and X are assumed to be aligned");
#endif
    if( singleL11AllGather )
        internal::ForwardSingle( L, X );
    else
        internal::ForwardMany( L, X );
}

template<typename F>
inline void FrontLowerForwardSolve( const DistMatrix<F>& L, DistMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerForwardSolve");
    if( L.Grid() != X.Grid() )
        LogicError("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
#endif
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const Int snSize = L.Width();
    DistMatrix<F> LT(g),
                  LB(g);
    LockedPartitionDown
    ( L, LT,
         LB, snSize );
    DistMatrix<F> XT(g),
                  XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    // XT := LT XT
    // TODO: Replace with TrsmLLNMedium?
    elem::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), LT, XT );

    // XB := XB - LB XT
    elem::Gemm( NORMAL, NORMAL, F(-1), LB, XT, F(1), XB );
}

template<typename F>
inline void FrontLowerBackwardSolve
( Orientation orientation, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool singleL11AllGather )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerBackwardSolve");
    if( L.Grid() != X.Grid() )
        LogicError("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( L.ColAlign() != X.ColAlign() )
        LogicError("L and X are assumed to be aligned");
    if( orientation == NORMAL )
        LogicError("This solve must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( orientation, L.LockedMatrix(), X.Matrix() );
        return;
    }

    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, L.Width() );

    if( XB.Height() != 0 )
    {
        // Subtract off the parent updates
        DistMatrix<F,STAR,STAR> Z(g);
        elem::LocalGemm( orientation, NORMAL, F(-1), LB, XB, Z );
        XT.SumScatterUpdate( F(1), Z );
    }

    if( singleL11AllGather )
        internal::BackwardSingle( orientation, LT, XT );
    else
        internal::BackwardMany( orientation, LT, XT );
}

template<typename F>
inline void FrontLowerBackwardSolve
( Orientation orientation, const DistMatrix<F>& L, DistMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerBackwardSolve");
    if( L.Grid() != X.Grid() )
        LogicError("L and X must be distributed over the same grid");
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
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        FrontLowerBackwardSolve( orientation, L.LockedMatrix(), X.Matrix() );
        return;
    }

    DistMatrix<F> LT(g), LB(g);
    LockedPartitionDown( L, LT, LB, L.Width() );

    DistMatrix<F> XT(g), XB(g);
    PartitionDown( X, XT, XB, L.Width() );

    elem::Gemm( orientation, NORMAL, F(-1), LB, XB, F(1), XT );
    elem::Trsm( LEFT, LOWER, orientation, NON_UNIT, F(1), LT, XT );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERSOLVE_DISTFRONT_HPP
