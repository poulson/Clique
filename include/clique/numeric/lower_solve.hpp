/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
    CallStackEntry entry("LowerSolve");
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
}

} // namespace cliq
