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
void clique::numeric::SupernodeLDLDiagonalSolve
( int supernodeSize,
  F alpha, const DistMatrix<F,VC,STAR>& d, 
                 DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::SupernodeLDLDiagonalSolve");
    if( d.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("Incompatible alignments");
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    if( d.Height() != X.Height() )
        throw std::logic_error("Invalid height of X");
#endif
    basic::DiagonalSolve
    ( LEFT, NORMAL, d.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::SupernodeLDLForwardSolve
( int supernodeSize,
  F alpha, const DistMatrix<F,VC,STAR>& L, 
                 DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::SupernodeLDLForwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
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
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
#endif
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
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    basic::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XT.Height() < supernodeSize )
    {
        const int blocksize = std::min(Blocksize(),supernodeSize-X0.Height());
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
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, (F)1, L11_STAR_STAR, X1_STAR_STAR );

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, L21, X1_STAR_STAR, (F)1, X2 );
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
void clique::numeric::SupernodeLDLBackwardSolve
( Orientation orientation, 
  int supernodeSize,
  F alpha, const DistMatrix<F,VC,STAR>& U,
                 DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::SupernodeLDLBackwardSolve");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( U.Height() != L.Width() || U.Height() != X.Height() || 
        U.Height() < supernodeSize )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  supernodeSize ~ " << supernodeSize << "\n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( U.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("U and X are assumed to be aligned");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    basic::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, U.Height()-supernodeSize );
    PartitionUp
    ( X, XT,
         XB, X.Height()-supernodeSize );

    // Subtract off the parent updates
    if( XT.Height() >= XB.Height() )
    {
        DistMatrix<F,STAR,STAR> XB_STAR_STAR(g);
        XB_STAR_STAR = XB;
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, UTR, XB_STAR_STAR, (F)1, XT );
    }
    else
    {
        DistMatrix<F,STAR,STAR> ZT(g);
        ZT.ResizeTo( XT.Height(), XT.Width() );
        basic::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, UBL, XB, (F)0, ZT );
        // Need SumScatter for [* ,* ] -> [VC,* ]
        // XT.SumScatterUpdate( (F)-1, ZT );
        throw std::logic_error("This update type is not yet finished");
    }

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
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := U11^-1[* ,* ] X1[* ,* ]
        basic::internal::LocalTrsm
        ( LEFT, UPPER, NORMAL, UNIT, (F)1, U11_STAR_STAR, X1_STAR_STAR );
        X1 = X1_STAR_STAR;

        // X0[VC,* ] -= U01[VC,* ] X1[* ,* ]
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, U01, X1_STAR_STAR, (F)1, X0 );
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

template void clique::numeric::SupernodeLDLForwardSolve
( int supernodeSize,
  float alpha, const DistMatrix<float,VC,STAR>& L,
                     DistMatrix<float,VC,STAR>& X );
template void clique::numeric::SupernodeLDLDiagonalSolve
( int supernodeSize,
  float alpha, const DistMatrix<float,VC,STAR>& d, 
                     DistMatrix<float,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::SupernodeLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  float alpha, const DistMatrix<float,VC,STAR>& U,
                     DistMatrix<float,VC,STAR>& X );

template void clique::numeric::SupernodeLDLForwardSolve
( int supernodeSize,
  double alpha, const DistMatrix<double,VC,STAR>& L, 
                      DistMatrix<double,VC,STAR>& X );
template void clique::numeric::SupernodeLDLDiagonalSolve
( int supernodeSize,
  double alpha, const DistMatrix<double,VC,STAR>& d, 
                      DistMatrix<double,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::SupernodeLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  double alpha, const DistMatrix<double,VC,STAR>& U,
                      DistMatrix<double,VC,STAR>& X );

template void clique::numeric::SupernodeLDLForwardSolve
( int supernodeSize,
  std::complex<float> alpha, const DistMatrix<std::complex<float>,VC,STAR>& L, 
                                   DistMatrix<std::complex<float>,VC,STAR>& X );
template void clique::numeric::SupernodeLDLDiagonalSolve
( int supernodeSize,
  std::complex<float> alpha, const DistMatrix<std::complex<float>,VC,STAR>& d, 
                                   DistMatrix<std::complex<float>,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::SupernodeLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  std::complex<float> alpha, const DistMatrix<std::complex<float>,VC,STAR>& U, 
                                   DistMatrix<std::complex<float>,VC,STAR>& X );

template void clique::numeric::SupernodeLDLForwardSolve
( int supernodeSize,
  std::complex<double> alpha, const DistMatrix<std::complex<double>,VC,STAR>& L, 
                                    DistMatrix<std::complex<double>,VC,STAR>& X );
template void clique::numeric::SupernodeLDLDiagonalSolve
( int supernodeSize,
  std::complex<double> alpha, const DistMatrix<std::complex<double>,VC,STAR>& d, 
                                    DistMatrix<std::complex<double>,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::SupernodeLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  std::complex<double> alpha, const DistMatrix<std::complex<double>,VC,STAR>& U,
                                    DistMatrix<std::complex<double>,VC,STAR>& X );
