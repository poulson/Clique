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
#ifndef CLIQ_NUMERIC_LDL_HPP
#define CLIQ_NUMERIC_LDL_HPP

namespace cliq {

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling the following factorizations.

template<typename F>
void InitializeDistLeaf( const DistSymmInfo& info, DistSymmFrontTree<F>& L );

template<typename F>
void LDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, 
  SymmFrontType newFrontType=LDL_2D );

} // namespace cliq

// Implementation
// ==============

#include "./ldl/local_front.hpp"
#include "./ldl/local_front_block.hpp"
#include "./ldl/dist_front.hpp"
#include "./ldl/dist_front_block.hpp"

#include "./ldl/local.hpp"
#include "./ldl/dist.hpp"

namespace cliq {

template<typename F>
inline void InitializeDistLeaf
( const DistSymmInfo& info, DistSymmFrontTree<F>& L )
{
    DEBUG_ONLY(CallStackEntry cse("InitializeDistLeaf"))
    const DistSymmNodeInfo& node = info.distNodes[0];
    Matrix<F>& topLocalFrontL = L.localFronts.back().frontL;
    DistMatrix<F>& front2dL = L.distFronts[0].front2dL;

    front2dL.LockedAttach
    ( topLocalFrontL.Height(), topLocalFrontL.Width(), *node.grid, 0, 0,
      topLocalFrontL.LockedBuffer(), topLocalFrontL.LDim() );
}

template<typename F>
inline void 
LDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    if( !Unfactored(L.frontType) )
        LogicError("Matrix is already factored");

    // Convert from 1D to 2D if necessary
    ChangeFrontType( L, SYMM_2D );

    // Perform the initial factorization
    L.frontType = InitialFactorType(newFrontType);
    LocalLDL( info, L );
    DistLDL( info, L );

    // Convert the fronts from the initial factorization to the requested form
    ChangeFrontType( L, newFrontType );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_HPP
