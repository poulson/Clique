/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F> 
void LocalFrontLDL( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalFrontLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("LocalFrontLDL");
    if( ABR.Height() != ABR.Width() )
        throw std::logic_error("ABR must be square");
    if( AL.Height() != AL.Width() + ABR.Width() )
        throw std::logic_error("AL and ABR don't have conformal dimensions");
    if( orientation == NORMAL )
        throw std::logic_error("LocalFrontLDL must be (conjugate-)transposed.");
#endif
    Matrix<F>
        ALTL, ALTR,  AL00, AL01, AL02,
        ALBL, ALBR,  AL10, AL11, AL12,
                     AL20, AL21, AL22;
    Matrix<F> d1;
    Matrix<F> S21;

    Matrix<F> S21T,
              S21B;
    Matrix<F> AL21T,
              AL21B;

    // Start the algorithm
    elem::PartitionDownDiagonal
    ( AL, ALTL, ALTR,
          ALBL, ALBR, 0 );
    while( ALTL.Width() < AL.Width() )
    {
        elem::RepartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, /**/ AL01, AL02,
         /***************/ /*********************/
                /**/        AL10, /**/ AL11, AL12,
          ALBL, /**/ ALBR,  AL20, /**/ AL21, AL22 );

        //--------------------------------------------------------------------//
        // This routine is unblocked, hence the need for us to generalize to 
        // an (ideally) faster blocked algorithm.
        elem::internal::LDLVar3( orientation, AL11, d1 );

        elem::Trsm( RIGHT, LOWER, orientation, UNIT, F(1), AL11, AL21 );

        S21 = AL21;
        elem::DiagonalSolve( RIGHT, NORMAL, d1, AL21 );

        elem::PartitionDown
        ( S21, S21T,
               S21B, AL22.Width() );
        elem::PartitionDown
        ( AL21, AL21T,
                AL21B, AL22.Width() );
        elem::Gemm( NORMAL, orientation, F(-1), S21, AL21T, F(1), AL22 );
        elem::internal::TrrkNT
        ( LOWER, orientation, F(-1), S21B, AL21B, F(1), ABR );
        //--------------------------------------------------------------------//

        elem::SlidePartitionDownDiagonal
        ( ALTL, /**/ ALTR,  AL00, AL01, /**/ AL02,
                /**/        AL10, AL11, /**/ AL12,
         /***************/ /*********************/
          ALBL, /**/ ALBR,  AL20, AL21, /**/ AL22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
