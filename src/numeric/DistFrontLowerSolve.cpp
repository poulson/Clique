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

template<typename F>
void clique::numeric::DistFrontLowerForwardSolve
( Diagonal diag, const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLowerForwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
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
        LocalFrontLowerForwardSolve
        ( diag, L.LockedLocalMatrix(), X.LocalMatrix(), checkIfSingular );
#ifndef RELEASE
        PopCallStack();
#endif
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

    elemental::LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    elemental::PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        elemental::LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        elemental::RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        elemental::basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, (F)1, L11_STAR_STAR, X1_STAR_STAR,
          checkIfSingular );
        X1 = X1_STAR_STAR;

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        elemental::basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, L21, X1_STAR_STAR, (F)1, X2 );
        //--------------------------------------------------------------------//

        elemental::SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        elemental::SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::DistFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::DistFrontLowerBackwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( orientation == NORMAL )
        throw std::logic_error("This solve must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();
    if( g.Size() == 1 )
    {
        LocalFrontLowerBackwardSolve
        ( orientation, diag, L.LockedLocalMatrix(), X.LocalMatrix(),
          checkIfSingular );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    elemental::LockedPartitionDown
    ( L, LT,
         LB, L.Width() );

    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    elemental::PartitionDown
    ( X, XT,
         XB, L.Width() );

    if( XB.Height() != 0 )
    {
        // Subtract off the parent updates
        DistMatrix<F,STAR,STAR> Z(XT.Height(),XT.Width(),g);
        Z.ResizeTo( XT.Height(), XT.Width() );
        elemental::basic::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, LB, XB, (F)0, Z );
        XT.SumScatterUpdate( (F)1, Z );
    }

    elemental::basic::internal::TrsmLLTSmall
    ( orientation, diag, (F)1, LT, XT, checkIfSingular );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::DistFrontLowerForwardSolve
( Diagonal diag, 
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X,
        bool checkIfSingular );
template void clique::numeric::DistFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X,
        bool checkIfSingular );

template void clique::numeric::DistFrontLowerForwardSolve
( Diagonal diag,
  const DistMatrix<double,VC,STAR>& L, 
        DistMatrix<double,VC,STAR>& X,
        bool checkIfSingular );
template void clique::numeric::DistFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const DistMatrix<double,VC,STAR>& L,
        DistMatrix<double,VC,STAR>& X,
        bool checkIfSingular );

template void clique::numeric::DistFrontLowerForwardSolve
( Diagonal diag,
  const DistMatrix<std::complex<float>,VC,STAR>& L, 
        DistMatrix<std::complex<float>,VC,STAR>& X,
        bool checkIfSingular );
template void clique::numeric::DistFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const DistMatrix<std::complex<float>,VC,STAR>& L, 
        DistMatrix<std::complex<float>,VC,STAR>& X,
        bool checkIfSingular );

template void clique::numeric::DistFrontLowerForwardSolve
( Diagonal diag,
  const DistMatrix<std::complex<double>,VC,STAR>& L, 
        DistMatrix<std::complex<double>,VC,STAR>& X,
        bool checkIfSingular );
template void clique::numeric::DistFrontLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  const DistMatrix<std::complex<double>,VC,STAR>& L,
        DistMatrix<std::complex<double>,VC,STAR>& X,
        bool checkIfSingular );
