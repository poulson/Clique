/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T>
void DistFrontLowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X );

template<typename T>
void DistFrontLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X );

template<typename T>
void DistFrontLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename T>
void ModifyForTrmm
( DistMatrix<T,STAR,STAR>& D, UnitOrNonUnit diag, int diagOffset )
{
#ifndef RELEASE
    cliq::PushCallStack("ModifyForTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        MemZero( D.LocalBuffer(j,j), length );
        if( diag == UNIT && j-diagOffset < height )
            D.SetLocal( j-diagOffset, j, T(1) );
    }
#ifndef RELEASE
    cliq::PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void DistFrontLowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("DistFrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        DistFrontLowerMultiplyNormal( diag, diagOffset, L, X );
    else
        DistFrontLowerMultiplyTranspose( orientation, diag, diagOffset, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void DistFrontLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("DistFrontLowerMultiplyNormal");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
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
    elem::LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, L.Width() );
    elem::PartitionDown
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
        elem::internal::LocalGemm
        ( NORMAL, NORMAL, T(1), L21, X1_STAR_STAR, T(1), X2 );

        if( diagOffset == 0 )
        {
            L11_STAR_STAR = L11;
            elem::internal::LocalTrmm
            ( LEFT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            L11_STAR_STAR = L11;
            internal::ModifyForTrmm( L11_STAR_STAR, diag, diagOffset );
            elem::internal::LocalTrmm
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void DistFrontLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistMatrix<T,VC,STAR>& L, DistMatrix<T,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("DistFrontLowerMultiplyTranspose");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transposed");
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
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
        if( diagOffset == 0 )
        {
            elem::internal::LocalTrmm
            ( LEFT, LOWER, orientation, diag, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        else
        {
            internal::ModifyForTrmm( L11_STAR_STAR, diag, diagOffset );
            elem::internal::LocalTrmm
            ( LEFT, LOWER, orientation, NON_UNIT, 
              T(1), L11_STAR_STAR, X1_STAR_STAR );
        }
        X1 = X1_STAR_STAR;

        Z1_STAR_STAR.ResizeTo( X1.Height(), X1.Width() );
        elem::internal::LocalGemm
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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
