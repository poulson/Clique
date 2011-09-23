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
void clique::numeric::SetSolveMode( DistSymmFact<F>& distL, SolveMode mode )
{
#ifndef RELEASE
    PushCallStack("numeric::SetSolveMode");
#endif
    // Check if this call can be a no-op
    if( mode == distL.mode ) 
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    distL.mode = mode;
    const int numSupernodes = distL.supernodes.size();    
    if( numSupernodes == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    DistSymmFactSupernode<F>& leafSN = distL.supernodes[0];
    if( mode == FEW_RHS )
    {
        leafSN.front1d.LockedView
        ( leafSN.front2d.Height(), leafSN.front2d.Width(), 0,
          leafSN.front2d.LockedLocalBuffer(), leafSN.front2d.LocalLDim(),
          leafSN.front2d.Grid() );
        for( int s=1; s<numSupernodes; ++s )
        {
            DistSymmFactSupernode<F>& sn = distL.supernodes[s];
            sn.front1d.Empty();
            sn.front1d.SetGrid( sn.front2d.Grid() );
            sn.front1d = sn.front2d;
            sn.front2d.Empty();
        }
    }
    else
    {
        leafSN.front2d.LockedView
        ( leafSN.front1d.Height(), leafSN.front1d.Width(), 0, 0,
          leafSN.front1d.LockedLocalBuffer(), leafSN.front1d.LocalLDim(),
          leafSN.front1d.Grid() );
        for( int s=1; s<numSupernodes; ++s )
        {
            DistSymmFactSupernode<F>& sn = distL.supernodes[s];
            sn.front2d.Empty();
            sn.front2d.SetGrid( sn.front1d.Grid() );
            sn.front2d = sn.front1d;
            sn.front1d.Empty();
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::SetSolveMode
( DistSymmFact<float>& distL, SolveMode mode );

template void clique::numeric::SetSolveMode
( DistSymmFact<double>& distL, SolveMode mode );

template void clique::numeric::SetSolveMode
( DistSymmFact<std::complex<float> >& distL, SolveMode mode );

template void clique::numeric::SetSolveMode
( DistSymmFact<std::complex<double> >& distL, SolveMode mode );
