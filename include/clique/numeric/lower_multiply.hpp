/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T>
void LowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX );

} // namespace cliq

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./lower_multiply/local_front.hpp"
#include "./lower_multiply/dist_front.hpp"

#include "./lower_multiply/local.hpp"
#include "./lower_multiply/dist.hpp"

namespace cliq {

template<typename T>
inline void LowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )
{
#ifndef RELEASE
    PushCallStack("LowerMultiply");
#endif
    if( orientation == NORMAL )
    {
        LocalLowerMultiplyNormal( diag, diagOffset, info, L, localX );
        DistLowerMultiplyNormal( diag, diagOffset, info, L, localX );
    }
    else
    {
        DistLowerMultiplyTranspose
        ( orientation, diag, diagOffset, info, L, localX );
        LocalLowerMultiplyTranspose
        ( orientation, diag, diagOffset, info, L, localX );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
