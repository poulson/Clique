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

namespace internal {
template<typename F>
void AddInLocalData
( const DistMatrix<F,VC,STAR>& X1, DistMatrix<F,STAR,STAR>& Z )
{
#ifndef RELEASE
    PushCallStack("internal::AddInLocalData");
#endif
    const int width = X1.Width();
    const int localHeight = X1.LocalHeight();
    const int stride = X1.Grid().Size();
    const int offset = X1.ColShift();
    for( int j=0; j<width; ++j )
    {
        F* ZColBuffer = Z.LocalBuffer(0,j);
        const F* X1ColBuffer = X1.LockedLocalBuffer(0,j);
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            ZColBuffer[offset+stride*iLocal] += X1ColBuffer[iLocal];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
}

template<typename F>
void clique::numeric::DistFrontLDLDiagonalSolve
( int supernodeSize, const DistMatrix<F,VC,STAR>& d, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLDLDiagonalSolve");
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
void clique::numeric::DistFrontLDLForwardSolve
( int supernodeSize, const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLDLForwardSolve");
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
    if( g.Size() == 1 )
    {
        LocalFrontLDLForwardSolve
        ( supernodeSize, L.LockedLocalMatrix(), X.LocalMatrix() );
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
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, UNIT, (F)1, L11_STAR_STAR, X1_STAR_STAR );
        X1 = X1_STAR_STAR;

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
void clique::numeric::DistFrontLDLBackwardSolve
( Orientation orientation, 
  int supernodeSize, const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLDLBackwardSolve");
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
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        LocalFrontLDLBackwardSolve
        ( orientation, supernodeSize, L.LockedLocalMatrix(), X.LocalMatrix() );
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
    DistMatrix<F,STAR,STAR> Z(g);

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, supernodeSize );
    PartitionDown
    ( X, XT,
         XB, supernodeSize );

    // Subtract off the parent updates
    Z.ResizeTo( XT.Height(), XT.Width() );
    basic::internal::LocalGemm( orientation, NORMAL, (F)-1, LBL, XB, (F)0, Z );
    XT.SumScatterUpdate( (F)1, Z );
    Z.Empty();

    // Solve the remaining triangular system
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,   L00, L01, /**/ L02,
               /**/        L10, L11, /**/ L12,
         /*************/  /******************/
          LBL, /**/ LBR,   L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        // X1 -= L21' X2
        Z.ResizeTo( X1.Height(), X1.Width() );
        basic::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, L21, X2, (F)0, Z );
        internal::AddInLocalData( X1, Z );
        Z.SumOverGrid();

        // X1 := L11^-1 X1
        L11_STAR_STAR = L11; 
        basic::internal::LocalTrsm
        ( LEFT, LOWER, orientation, UNIT, (F)1, L11_STAR_STAR, Z );
        X1 = Z;
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
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::DistFrontLDLForwardSolve
( int supernodeSize,
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X );
template void clique::numeric::DistFrontLDLDiagonalSolve
( int supernodeSize,
  const DistMatrix<float,VC,STAR>& d, 
        DistMatrix<float,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::DistFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X );

template void clique::numeric::DistFrontLDLForwardSolve
( int supernodeSize,
  const DistMatrix<double,VC,STAR>& L, 
        DistMatrix<double,VC,STAR>& X );
template void clique::numeric::DistFrontLDLDiagonalSolve
( int supernodeSize,
  const DistMatrix<double,VC,STAR>& d, 
        DistMatrix<double,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::DistFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const DistMatrix<double,VC,STAR>& L,
        DistMatrix<double,VC,STAR>& X );

template void clique::numeric::DistFrontLDLForwardSolve
( int supernodeSize,
  const DistMatrix<std::complex<float>,VC,STAR>& L, 
        DistMatrix<std::complex<float>,VC,STAR>& X );
template void clique::numeric::DistFrontLDLDiagonalSolve
( int supernodeSize,
  const DistMatrix<std::complex<float>,VC,STAR>& d, 
        DistMatrix<std::complex<float>,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::DistFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const DistMatrix<std::complex<float>,VC,STAR>& L, 
        DistMatrix<std::complex<float>,VC,STAR>& X );

template void clique::numeric::DistFrontLDLForwardSolve
( int supernodeSize,
  const DistMatrix<std::complex<double>,VC,STAR>& L, 
        DistMatrix<std::complex<double>,VC,STAR>& X );
template void clique::numeric::DistFrontLDLDiagonalSolve
( int supernodeSize,
  const DistMatrix<std::complex<double>,VC,STAR>& d, 
        DistMatrix<std::complex<double>,VC,STAR>& X,
  bool checkIfSingular );
template void clique::numeric::DistFrontLDLBackwardSolve
( Orientation orientation, int supernodeSize,
  const DistMatrix<std::complex<double>,VC,STAR>& L,
        DistMatrix<std::complex<double>,VC,STAR>& X );
