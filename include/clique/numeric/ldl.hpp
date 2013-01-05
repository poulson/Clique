/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
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

    front2dL.LockedAttach
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
