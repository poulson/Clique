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
void clique::numeric::LocalLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
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
        WB.SetToZero();

        // Update using the children (if they exist)
        const int numChildren = sn.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = sn.children[0];
            const int rightIndex = sn.children[1];
            Matrix<F>& leftWork = L.local.fronts[leftIndex].work;
            Matrix<F>& rightWork = L.local.fronts[rightIndex].work;
            const int leftSupernodeSize = S.local.supernodes[leftIndex].size;
            const int rightSupernodeSize = S.local.supernodes[rightIndex].size;
            const int leftUpdateSize = leftWork.Height()-leftSupernodeSize;
            const int rightUpdateSize = rightWork.Height()-rightSupernodeSize;

            // Add the left child's update onto ours
            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftWork, leftSupernodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = sn.leftChildRelIndices[iChild]; 
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightWork, rightSupernodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = sn.rightChildRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Solve against this front
        W.Print("W before forward solve");
        LocalFrontLDLForwardSolve( sn.size, front, W );
        W.Print("W after forward solve");

        // Store the supernode portion of the result
        XT = WT;
        XT.Print("XT after forward solve");
    }

    // Ensure that all of the temporary buffers are freed (except the root)
    for( int s=0; s<numSupernodes-1; ++s )
        L.local.fronts[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X, bool checkIfSingular )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLDiagonalSolve");
#endif
    const int numSupernodes = S.local.supernodes.size();
    const int width = X.Width();
    Matrix<F> XSub;
    for( int s=0; s<numSupernodes; ++s )
    {
        const LocalSymmFactSupernode& sn = S.local.supernodes[s];
        const Matrix<F>& front = L.local.fronts[s].front;
        XSub.View( X, sn.myOffset, 0, sn.size, width );

        Matrix<F> frontTL;
        frontTL.LockedView( front, 0, 0, sn.size, sn.size );
        Matrix<F> d;
        frontTL.GetDiagonal( d );
        XSub.Print("XSub before diagonal solve");
        LocalFrontLDLDiagonalSolve( sn.size, d, XSub, checkIfSingular );
        XSub.Print("XSub after diagonal solve");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S, 
  const numeric::SymmFrontTree<F>& L,
        Matrix<F>& X )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLBackwardSolve");
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

        // Update using the parent
        const int parent = sn.parent;
        Matrix<F>& parentWork = L.local.fronts[parent].work;
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

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( sn.isLeftChild )
            parentWork.Empty();

        // Solve against this front
        W.Print("W before backward solve");
        LocalFrontLDLBackwardSolve( orientation, sn.size, front, W );
        W.Print("W after backward solve");

        // Store the supernode portion of the result
        XT = WT;
        XT.Print("XT after backward solve");
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    for( int s=0; s<numSupernodes; ++s )
        L.local.fronts[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<float>& L,
        Matrix<float>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<double>& L,
        Matrix<double>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<float> >& L,
        Matrix<std::complex<float> >& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::SymmFact& S,
  const numeric::SymmFrontTree<std::complex<double> >& L,
        Matrix<std::complex<double> >& X );
