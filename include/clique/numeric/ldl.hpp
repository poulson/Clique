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
#ifndef CLIQUE_LDL_HPP
#define CLIQUE_LDL_HPP 1

namespace cliq {

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling the following factorizations.

template<typename F>
void InitializeDistLeaf( const SymmInfo& info, SymmFrontTree<F>& L );

template<typename F>
void LDL( Orientation orientation, SymmInfo& info, SymmFrontTree<F>& L );

} // namespace cliq

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./ldl/local_front.hpp"
#include "./ldl/dist_front.hpp"

#include "./ldl/local.hpp"
#include "./ldl/dist.hpp"

namespace cliq {

template<typename F>
inline void InitializeDistLeaf( const SymmInfo& info, SymmFrontTree<F>& L )
{
#ifndef RELEASE
    PushCallStack("InitializeDistLeaf");
#endif
    const DistSymmNodeInfo& node = info.dist.nodes[0];
    Matrix<F>& topLocalFrontL = L.local.fronts.back().frontL;
    DistMatrix<F>& front2dL = L.dist.fronts[0].front2dL;

    front2dL.LockedView
    ( topLocalFrontL.Height(), topLocalFrontL.Width(), 0, 0,
      topLocalFrontL.LockedBuffer(), topLocalFrontL.LDim(), 
      *node.grid );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void LDL( Orientation orientation, SymmInfo& info, SymmFrontTree<F>& L )
{
#ifndef RELEASE
    PushCallStack("LDL");
#endif
    LocalLDL( orientation, info, L );
    DistLDL( orientation, info, L );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif /* CLIQUE_LDL_HPP */
