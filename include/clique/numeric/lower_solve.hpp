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
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F>
void LowerSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

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
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("LowerSolve");
#endif
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( info, L, X );
        DistLowerForwardSolve( info, L, X );
    }
    else
    {
        DistLowerBackwardSolve( orientation, info, L, X );
        LocalLowerBackwardSolve( orientation, info, L, X );
    }
}

template<typename F>
inline void LowerSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("LowerSolve");
#endif
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( info, L, X );
        DistLowerForwardSolve( info, L, X );
    }
    else
    {
        DistLowerBackwardSolve( orientation, info, L, X );
        LocalLowerBackwardSolve( orientation, info, L, X );
    }
}

} // namespace cliq
