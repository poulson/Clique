/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F>
void LowerSolve
( Orientation orientation, UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

} // namespace cliq

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./lower_solve/local_front.hpp"
#include "./lower_solve/local_front_block.hpp"
#include "./lower_solve/dist_front.hpp"
#include "./lower_solve/dist_front_fast.hpp"
#include "./lower_solve/dist_front_block.hpp"

#include "./lower_solve/local.hpp"
#include "./lower_solve/dist.hpp"

namespace cliq {

template<typename F>
inline void LowerSolve
( Orientation orientation, UnitOrNonUnit diag, 
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("LowerSolve");
#endif
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( diag, info, L, localX );
        DistLowerForwardSolve( diag, info, L, localX );
    }
    else
    {
        DistLowerBackwardSolve( orientation, diag, info, L, localX );
        LocalLowerBackwardSolve( orientation, diag, info, L, localX );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
