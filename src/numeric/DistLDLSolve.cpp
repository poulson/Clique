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

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<F>& localL,
  const numeric::DistSymmFact<F>& distL,
        Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
#endif
    const int numSupernodes = distS.supernodes.size();
    const int width = localX.Width();
    if( distL.mode == MANY_RHS )
        throw std::logic_error("This solve mode is not yet implemented");
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Copy the information from the local portion into the distributed root
    const LocalSymmFactSupernode<F>& topLocalSN = localL.supernodes.back();
    const DistSymmFactSupernode<F>& bottomDistSN = distL.supernodes[0];
    bottomDistSN.workspace2d.LocalMatrix().LockedView( topLocalSN.workspace );
    
    // Perform the distributed portion of the forward solve
    for( int k=1; k<numSupernodes; ++k )
    {
        const symbolic::DistSymmFactSupernode& symbSN = distS.supernodes[k];
        const numeric::DistSymmFactSupernode<F>& numSN = distL.supernodes[k];
        const Grid& grid = numSN.front1d.Grid();

        // Set up a workspace
        DistMatrix<F,VC,STAR>& W = numSN.workspace1d;
        W.ResizeTo( numSN.front1d.Height(), width );
        DistMatrix<F,VC,STAR> WT(grid), WB(grid);
        WT.View( W, 0, 0, symbSN.size, width );
        WB.View( W, symbSN.size, 0, W.Height()-symbSN.size, width );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.LockedView
        ( localX, symbSN.localOffset1d, 0, symbSN.localSize1d, width );
        WT.LocalMatrix() = XT;
        WB.SetToZero();

        // Update using the children 
        // HERE
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<float>& localL,
  const numeric::DistSymmFact<float>& distL,
        Matrix<float>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<double>& localL,
  const numeric::DistSymmFact<double>& distL,
        Matrix<double>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
  const numeric::DistSymmFact<std::complex<float> >& distL,
        Matrix<std::complex<float> >& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
  const numeric::DistSymmFact<std::complex<double> >& distL,
        Matrix<std::complex<double> >& localX );
