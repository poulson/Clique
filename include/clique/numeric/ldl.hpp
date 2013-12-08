/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
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

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

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
    ( topLocalFrontL.Height(), topLocalFrontL.Width(), 0, 0,
      topLocalFrontL.LockedBuffer(), topLocalFrontL.LDim(), 
      *node.grid );
}

template<typename F>
inline void 
LDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    if( L.frontType != SYMM_2D )
        LogicError
        ("Should only perform LDL factorization of 2D symmetric/Hermitian "
         "matrices");

    if( newFrontType == LDL_2D        || newFrontType == LDL_1D || 
        newFrontType == LDL_SELINV_2D || newFrontType == LDL_SELINV_1D )
    {
        L.frontType = LDL_2D;
    }
    else if( newFrontType == LDL_INTRAPIV_2D )
    {
        L.frontType = LDL_INTRAPIV_2D;    
    }
    else if( newFrontType == BLOCK_LDL_2D || 
             newFrontType == BLOCK_LDL_INTRAPIV_2D )
    {
        L.frontType = newFrontType;
    }
    else
        LogicError("Invalid new front type");

    LocalLDL( info, L );
    DistLDL( info, L );

    ChangeFrontType( L, newFrontType );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_HPP
