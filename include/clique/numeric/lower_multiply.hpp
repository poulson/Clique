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
#ifndef CLIQUE_LOWER_MULTIPLY_HPP
#define CLIQUE_LOWER_MULTIPLY_HPP 1

namespace cliq {

template<typename F>
void LowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX );

} // namespace cliq

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./lower_multiply/local_front.hpp"
#include "./lower_multiply/dist_front.hpp"

#include "./lower_multiply/local.hpp"
#include "./lower_multiply/dist.hpp"

namespace cliq {

template<typename F>
inline void LowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("LowerMultiply");
#endif
    if( orientation == NORMAL )
    {
        LocalLowerMultiplyNormal( diag, diagOffset, info, L, localX );
        DistLowerMultiplyNormal( diag, diagOffset, info, L, localX );
    }
    else
    {
        DistLowerMultiplyTranspose
        ( orientation, diag, diagOffset, info, L, localX );
        LocalLowerMultiplyTranspose
        ( orientation, diag, diagOffset, info, L, localX );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // CLIQUE_LOWER_MULTIPLY_HPP
