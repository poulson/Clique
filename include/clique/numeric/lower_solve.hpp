/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERSOLVE_HPP
#define CLIQ_NUMERIC_LOWERSOLVE_HPP

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
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( info, L, X );
        DistLowerForwardSolve( info, L, X );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DistLowerBackwardSolve( info, L, X, conjugate );
        LocalLowerBackwardSolve( info, L, X, conjugate );
    }
}

template<typename F>
inline void LowerSolve
( Orientation orientation, const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LowerSolve"))
    if( orientation == NORMAL )
    {
        LocalLowerForwardSolve( info, L, X );
        DistLowerForwardSolve( info, L, X );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DistLowerBackwardSolve( info, L, X, conjugate );
        LocalLowerBackwardSolve( info, L, X, conjugate );
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERSOLVE_HPP
