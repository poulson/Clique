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
void clique::numeric::LocalLowerMultiplyNormal
( Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalLowerMultiplyNormal");
#endif
    const int numSupernodes = S.local.supernodes.size();
    const int width = X.Width();
    for( int s=0; s<numSupernodes; ++s )
    {
        const LocalSymmFactSupernode& sn = S.local.supernodes[s];
        const Matrix<F>& front = L.local.fronts[s].front;
        Matrix<F>& W = L.local.fronts[s].work;

        // Set up a workspace
        W.ResizeTo( front.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, sn.myOffset, 0, sn.size, width );
        WT = XT;

        // Multiply with this front
        //LocalFrontLowerMultiply
        //( NORMAL, diag, diagOffset, sn.size, front, W );

        // Add the updates from the children
        throw std::logic_error("This routine is not finished");

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed (except the root)
    for( int s=0; s<numSupernodes-1; ++s )
        L.local.fronts[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S, 
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalLowerMultiplyTranspose");
#endif
    const int numSupernodes = S.local.supernodes.size();
    const int width = X.Width();
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Pull in the top local information from the bottom distributed information
    L.local.fronts.back().work.LockedView
    ( L.dist.fronts[0].work1d.LocalMatrix() );

    for( int s=numSupernodes-2; s>=0; --s )
    {
        const LocalSymmFactSupernode& sn = S.local.supernodes[s];
        const Matrix<F>& front = L.local.fronts[s].front;
        Matrix<F>& W = L.local.fronts[s].work;

        // Set up a workspace
        W.ResizeTo( front.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown
        ( W, WT,
             WB, sn.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, sn.myOffset, 0, sn.size, width );
        WT = XT;

        // Update using the parent's portion of the RHS
        const int parent = sn.parent;
        Matrix<F>& parentWork = L.local.fronts[parent].work;
        /*
        const LocalSymmFactSupernode& parentSN = S.local.supernodes[parent];
        const int currentUpdateSize = WB.Height();

        const std::vector<int>& parentRelIndices = 
          ( sn.isLeftChild ? 
            parentSN.leftChildRelIndices :
            parentSN.rightChildRelIndices );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelIndices[iCurrent];
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }
        */
        throw std::logic_error("This routine is not finished");

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( sn.isLeftChild )
            parentWork.Empty();

        // Multiply against this front
        //LocalFrontLowerMultiply
        //( orientation, diag, diagOffset, sn.size, front, W );

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    for( int s=0; s<numSupernodes; ++s )
        L.local.fronts[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLowerMultiplyNormal
( Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& X );
template void clique::numeric::LocalLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& X );

template void clique::numeric::LocalLowerMultiplyNormal
( Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& X );
template void clique::numeric::LocalLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& X );

template void clique::numeric::LocalLowerMultiplyNormal
( Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& X );
template void clique::numeric::LocalLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& X );

template void clique::numeric::LocalLowerMultiplyNormal
( Diagonal diag, int diagOffset, 
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& X );
template void clique::numeric::LocalLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& X );
