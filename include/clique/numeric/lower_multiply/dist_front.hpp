/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERMULTIPLY_DISTFRONT_HPP
#define CLIQ_NUMERIC_LOWERMULTIPLY_DISTFRONT_HPP

namespace cliq {

template<typename T>
void FrontLowerMultiply
( Orientation orientation, 
  int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X );

template<typename T>
void FrontLowerMultiplyNormal
( int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X );

template<typename T>
void FrontLowerMultiplyTranspose
( int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X,
  bool conjugate=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename T>
void ModifyForTrmm( DistMatrix<T,STAR,STAR>& D, int diagOff )
{
#ifndef RELEASE
    cliq::CallStackEntry cse("ModifyForTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOff,height-j);
        MemZero( D.Buffer(j,j), length );
    }
}

} // namespace internal

template<typename T>
inline void FrontLowerMultiply
( Orientation orientation, int diagOff,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        FrontLowerMultiplyNormal( diagOff, L, X );
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        FrontLowerMultiplyTranspose( diagOff, L, X, conjugate );
    }
}

template<typename T>
inline void FrontLowerMultiplyNormal
( int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerMultiplyNormal");
    if( L.Grid() != X.Grid() )
        LogicError("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( L.ColAlign() != X.ColAlign() )
        LogicError("L and X are assumed to be aligned");
    if( diagOff > 0 )
        LogicError("Diagonal offsets cannot be positive");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, L.Width() );
    PartitionDown
    ( X, XT,
         XB, L.Width() );
    while( XT.Height() > 0 )
    {
        elem::LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        elem::RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        X1_STAR_STAR = X1;
        elem::LocalGemm( NORMAL, NORMAL, T(1), L21, X1_STAR_STAR, T(1), X2 );

        if( diagOff == 0 )
        {
            L11_STAR_STAR = L11;
            elem::LocalTrmm
            ( LEFT, LOWER, NORMAL, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            L11_STAR_STAR = L11;
            internal::ModifyForTrmm( L11_STAR_STAR, diagOff );
            elem::LocalTrmm
            ( LEFT, LOWER, NORMAL, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;
        //--------------------------------------------------------------------//

        elem::SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        elem::SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
}

template<typename T>
inline void FrontLowerMultiplyTranspose
( int diagOff, const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X,
  bool conjugate )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerMultiplyTranspose");
    if( L.Grid() != X.Grid() )
        LogicError("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( L.ColAlign() != X.ColAlign() )
        LogicError("L and X are assumed to be aligned");
    if( diagOff > 0 )
        LogicError("Diagonal offsets cannot be positive");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,VC,STAR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,VC,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> X1_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> Z1_STAR_STAR(g);

    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( LTL.Width() < L.Width() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,   L00, /**/ L01, L02,
         /*************/  /******************/
               /**/        L10, /**/ L11, L12,
          LBL, /**/ LBR,   L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, L11.Height() );

        //--------------------------------------------------------------------//
        X1_STAR_STAR = X1; // Can this be avoided?
        L11_STAR_STAR = L11;
        if( diagOff == 0 )
        {
            elem::LocalTrmm
            ( LEFT, LOWER, orientation, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            internal::ModifyForTrmm( L11_STAR_STAR, diagOff );
            elem::LocalTrmm
            ( LEFT, LOWER, orientation, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;

        Z1_STAR_STAR.ResizeTo( X1.Height(), X1.Width() );
        elem::LocalGemm
        ( orientation, NORMAL, T(1), L21, X2, T(0), Z1_STAR_STAR );
        X1.SumScatterUpdate( T(1), Z1_STAR_STAR );
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

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERMULTIPLY_DISTFRONT_HPP
