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
#ifndef CLIQUE_BLOCK_LDL_HPP
#define CLIQUE_BLOCK_LDL_HPP 1

namespace cliq {

// All fronts of L are required to be initialized to the expansions of the 
// original sparse matrix before calling the following factorizations.

template<typename F>
void BlockLDL( Orientation orientation, SymmInfo& info, SymmFrontTree<F>& L );

} // namespace cliq

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./block_ldl/local_front.hpp"
#include "./block_ldl/dist_front.hpp"

#include "./block_ldl/local.hpp"
#include "./block_ldl/dist.hpp"

namespace cliq {

template<typename F>
inline void BlockLDL
( Orientation orientation, SymmInfo& info, SymmFrontTree<F>& L )
{
#ifndef RELEASE
    PushCallStack("BlockLDL");
#endif
    LocalBlockLDL( orientation, info, L );
    DistBlockLDL( orientation, info, L );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif /* CLIQUE_BLOCK_LDL_HPP */
