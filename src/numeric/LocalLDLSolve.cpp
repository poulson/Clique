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
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<F>& L,
        Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = X.Width();
    for( int s=0; s<numSupernodes; ++s )
    {
        const symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[s];
        const numeric::LocalSymmFactSupernode<F>& numSN = L.supernodes[s];

        // Set up a workspace
        Matrix<F>& W = numSN.work;
        W.ResizeTo( numSN.front.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown
        ( W, WT,
             WB, symbSN.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, symbSN.myOffset, 0, symbSN.size, width );
        WT = XT;
        WB.SetToZero();

        // Update using the children (if they exist)
        const int numChildren = symbSN.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = symbSN.children[0];
            const int rightIndex = symbSN.children[1];
            Matrix<F>& leftWork = L.supernodes[leftIndex].work;
            Matrix<F>& rightWork = L.supernodes[rightIndex].work;
            const int leftSupernodeSize = S.supernodes[leftIndex].size;
            const int rightSupernodeSize = S.supernodes[rightIndex].size;
            const int leftUpdateSize = leftWork.Height()-leftSupernodeSize;
            const int rightUpdateSize = rightWork.Height()-rightSupernodeSize;

            // Add the left child's update onto ours
            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftWork, leftSupernodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = symbSN.leftChildRelIndices[iChild]; 
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, -leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightWork, rightSupernodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = symbSN.rightChildRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, -rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Call the custom supernode forward solve
        LocalSupernodeLDLForwardSolve( symbSN.size, numSN.front, W );

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed (except the root)
    for( int s=0; s<numSupernodes-1; ++s )
        L.supernodes[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<F>& L,
  Matrix<F>& X, bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLDiagonalSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = X.Width();
    Matrix<F> XSub;
    for( int s=0; s<numSupernodes; ++s )
    {
        const symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[s];
        const LocalSymmFactSupernode<F>& numSN = L.supernodes[s];
        XSub.View( X, symbSN.myOffset, 0, symbSN.size, width );

        Matrix<F> frontTL;
        frontTL.LockedView( numSN.front, 0, 0, symbSN.size, symbSN.size );
        Matrix<F> d;
        frontTL.GetDiagonal( d );
        LocalSupernodeLDLDiagonalSolve( symbSN.size, d, XSub, checkIfSingular );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S, 
  const numeric::LocalSymmFact<F>& localL,
  const numeric::DistSymmFact<F>& distL,
        Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLBackwardSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = X.Width();
    if( numSupernodes == 0 || width == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Pull in the top local information from the bottom distributed information
    const LocalSymmFactSupernode<F>& topLocalSN = localL.supernodes.back();
    const DistSymmFactSupernode<F>& bottomDistSN = distL.supernodes[0];
    topLocalSN.work.LockedView( bottomDistSN.work1d.LocalMatrix() );

    for( int s=numSupernodes-2; s>=0; --s )
    {
        const symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[s];
        const numeric::LocalSymmFactSupernode<F>& numSN = localL.supernodes[s];

        // Set up a workspace
        Matrix<F>& W = numSN.work;
        W.ResizeTo( numSN.front.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown
        ( W, WT,
             WB, symbSN.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, symbSN.myOffset, 0, symbSN.size, width );
        WT = XT;

        // Update using the parent
        const int parentIndex = symbSN.parent;
        Matrix<F>& parentWork = localL.supernodes[parentIndex].work;
        const symbolic::LocalSymmFactSupernode& parentSymbSN = 
            S.supernodes[parentIndex];
        const int currentUpdateSize = WB.Height();

        const std::vector<int>& parentRelIndices = 
          ( symbSN.isLeftChild ? 
            parentSymbSN.leftChildRelIndices :
            parentSymbSN.rightChildRelIndices );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelIndices[iCurrent];
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( symbSN.isLeftChild )
            parentWork.Empty();

        // Call the custom supernode backward solve
        LocalSupernodeLDLBackwardSolve
        ( orientation, symbSN.size, numSN.front, W );

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    for( int s=0; s<numSupernodes; ++s )
        localL.supernodes[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<float>& L,
        Matrix<float>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<float>& L,
  Matrix<float>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<float>& localL,
  const numeric::DistSymmFact<float>& distL,
        Matrix<float>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<double>& L,
        Matrix<double>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<double>& L,
  Matrix<double>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<double>& localL,
  const numeric::DistSymmFact<double>& distL,
        Matrix<double>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& L,
        Matrix<std::complex<float> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& L,
  Matrix<std::complex<float> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
  const numeric::DistSymmFact<std::complex<float> >& distL,
        Matrix<std::complex<float> >& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& L,
        Matrix<std::complex<double> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& L,
        Matrix<std::complex<double> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
  const numeric::DistSymmFact<std::complex<double> >& distL,
        Matrix<std::complex<double> >& X );
