/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F> 
void LocalFrontBlockLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalFrontBlockLDL
( Orientation orientation, Matrix<F>& AL, Matrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("LocalFrontBlockLDL");
#endif
    Matrix<F> ATL,
              ABL;
    elem::PartitionDown
    ( AL, ATL,
          ABL, AL.Width() );
    
    // Make a copy of the original contents of ABL
    Matrix<F> BBL( ABL );

    // Call the standard routine
    LocalFrontLDL( orientation, AL, ABR );

    // Copy the original contents of ABL back
    ABL = BBL;

    // Overwrite ATL with inv(L D L^[T/H]) = L^[-T/H] D^{-1} L^{-1}
    elem::TriangularInverse( LOWER, UNIT, ATL );
    elem::Trdtrmm( orientation, LOWER, ATL );
    elem::MakeTrapezoidal( LEFT, LOWER, 0, ATL );
    if( orientation == TRANSPOSE )
    {
        Matrix<F> ATLTrans;
        elem::Transpose( ATL, ATLTrans );
        elem::MakeTrapezoidal( LEFT, UPPER, 1, ATLTrans );
        elem::Axpy( F(1), ATLTrans, ATL );
    }
    else
    {
        Matrix<F> ATLAdj;
        elem::Adjoint( ATL, ATLAdj );
        elem::MakeTrapezoidal( LEFT, UPPER, 1, ATLAdj );
        elem::Axpy( F(1), ATLAdj, ATL );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
