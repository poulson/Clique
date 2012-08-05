/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
#ifndef RELEASE
    PushCallStack("InitializeDistLeaf");
#endif
    const DistSymmNodeInfo& node = info.distNodes[0];
    Matrix<F>& topLocalFrontL = L.localFronts.back().frontL;
    DistMatrix<F>& front2dL = L.distFronts[0].front2dL;

    front2dL.LockedView
    ( topLocalFrontL.Height(), topLocalFrontL.Width(), 0, 0,
      topLocalFrontL.LockedBuffer(), topLocalFrontL.LDim(), 
      *node.grid );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
LDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType )
{
#ifndef RELEASE
    PushCallStack("LDL");
#endif
    if( L.frontType != SYMM_2D )
        throw std::logic_error
        ("Should only perform LDL factorization of 2D "
         "symmetric/Hermitian matrices");

    bool blockLDL;
    if( newFrontType == LDL_2D || newFrontType == LDL_1D || 
        newFrontType == LDL_SELINV_2D || newFrontType == LDL_SELINV_1D )
    {
        blockLDL = false;
        L.frontType = LDL_2D;
    }
    else if( newFrontType == BLOCK_LDL_2D )
    {
        blockLDL = true;
        L.frontType = BLOCK_LDL_2D;
    }
    else
        throw std::logic_error("Invalid new front type");

    LocalLDL( info, L, blockLDL );
    DistLDL( info, L, blockLDL );

    ChangeFrontType( L, newFrontType );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
