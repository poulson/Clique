/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
#include "clique.hpp"
using namespace elemental;

// This routine could be modified later so that it uses much less memory
// by replacing the '=' redistributions with piece-by-piece redistributions.
template<typename F>
void clique::numeric::SetSolveMode( SymmFrontTree<F>& L, SolveMode mode )
{
#ifndef RELEASE
    PushCallStack("numeric::SetSolveMode");
#endif
    // Check if this call can be a no-op
    if( mode == L.dist.mode ) 
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    L.dist.mode = mode;
    const int numSupernodes = L.dist.fronts.size();    
    if( numSupernodes == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    DistSymmFront<F>& leafFront = L.dist.fronts[0];
    if( mode == FEW_RHS )
    {
        leafFront.front1d.LockedView
        ( leafFront.front2d.Height(), leafFront.front2d.Width(), 0,
          leafFront.front2d.LockedLocalBuffer(), leafFront.front2d.LocalLDim(),
          leafFront.front2d.Grid() );
        for( int s=1; s<numSupernodes; ++s )
        {
            DistSymmFront<F>& front = L.dist.fronts[s];
            front.front1d.Empty();
            front.front1d.SetGrid( front.front2d.Grid() );
            front.front1d = front.front2d;
            front.front2d.Empty();
        }
    }
    else
    {
        leafFront.front2d.LockedView
        ( leafFront.front1d.Height(), leafFront.front1d.Width(), 0, 0,
          leafFront.front1d.LockedLocalBuffer(), leafFront.front1d.LocalLDim(),
          leafFront.front1d.Grid() );
        for( int s=1; s<numSupernodes; ++s )
        {
            DistSymmFront<F>& front = L.dist.fronts[s];
            front.front2d.Empty();
            front.front2d.SetGrid( front.front1d.Grid() );
            front.front2d = front.front1d;
            front.front1d.Empty();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::SetSolveMode
( SymmFrontTree<float>& L, SolveMode mode );

template void clique::numeric::SetSolveMode
( SymmFrontTree<double>& L, SolveMode mode );

template void clique::numeric::SetSolveMode
( SymmFrontTree<std::complex<float> >& L, SolveMode mode );

template void clique::numeric::SetSolveMode
( SymmFrontTree<std::complex<double> >& L, SolveMode mode );
