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

template<typename F>
void BlockLowerSolve
( Orientation orientation, 
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

} // namespace cliq

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./block_lower_solve/local_front.hpp"
#include "./block_lower_solve/dist_front.hpp"

#include "./block_lower_solve/local.hpp"
#include "./block_lower_solve/dist.hpp"

namespace cliq {

template<typename F>
inline void BlockLowerSolve
( Orientation orientation,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("BlockLowerSolve");
#endif
    if( orientation == NORMAL )
    {
        LocalBlockLowerForwardSolve( info, L, localX );
        DistBlockLowerForwardSolve( info, L, localX );
    }
    else
    {
        DistBlockLowerBackwardSolve( orientation, info, L, localX );
        LocalBlockLowerBackwardSolve( orientation, info, L, localX );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
