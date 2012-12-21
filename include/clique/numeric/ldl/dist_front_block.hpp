/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F> 
void DistFrontBlockLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void DistFrontBlockLDL
( Orientation orientation, DistMatrix<F>& AL, DistMatrix<F>& ABR )
{
#ifndef RELEASE
    PushCallStack("internal::DistFrontBlockLDL");
#endif
    const Grid& g = AL.Grid();
    DistMatrix<F> ATL(g),
                  ABL(g);
    elem::PartitionDown
    ( AL, ATL,
          ABL, AL.Width() );

    // Make a copy of the original contents of ABL
    DistMatrix<F> BBL( ABL );

    // Call the standard routine
    DistFrontLDL( orientation, AL, ABR );

    // Copy the original contents of ABL back
    ABL = BBL;

    // Overwrite ATL with inv(L D L^[T/H]) = L^[-T/H] D^{-1} L^{-1}
    elem::TriangularInverse( LOWER, UNIT, ATL );
    elem::Trdtrmm( orientation, LOWER, ATL );
    elem::MakeTrapezoidal( LEFT, LOWER, 0, ATL );
    if( orientation == TRANSPOSE )
    {
        DistMatrix<F> ATLTrans( g );
        elem::Transpose( ATL, ATLTrans );
        elem::MakeTrapezoidal( LEFT, UPPER, 1, ATLTrans );
        elem::Axpy( F(1), ATLTrans, ATL );
    }
    else
    {
        DistMatrix<F> ATLAdj( g );
        elem::Adjoint( ATL, ATLAdj );
        elem::MakeTrapezoidal( LEFT, UPPER, 1, ATLAdj );
        elem::Axpy( F(1), ATLAdj, ATL );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
