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
  F alpha, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = X.Width();
    for( int k=0; k<numSupernodes; ++k )
    {
        const symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[k];
        const numeric::LocalSymmFactSupernode<F>& numSN = L.supernodes[k];

        // Set up a workspace
        Matrix<F>& W = numSN.workspace;
        W.ResizeTo( numSN.front.Height(), width );
        Matrix<F> WT, WB;
        WT.View( W, 0, 0, symbSN.size, width );
        WB.View( W, symbSN.size, 0, W.Height()-symbSN.size, width );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.LockedView( X, symbSN.offset, 0, symbSN.size, X.Width() );
        WT = XT;
        WB.SetToZero();

        // Update using the children (if they exist)
        const int numChildren = symbSN.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = symbSN.children[0];
            const int rightIndex = symbSN.children[1];
            Matrix<F>& leftWork = L.supernodes[leftIndex].workspace;
            Matrix<F>& rightWork = L.supernodes[rightIndex].workspace;
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
                    W.Update( iFront, j, -leftUpdate.Get(iFront,j) );
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
                    W.Update( iFront, j, -rightUpdate.Get(iFront,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Call the custom supernode forward solve
        LocalSupernodeLDLForwardSolve( symbSN.size, alpha, numSN.front, W );

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed
    for( int k=0; k<numSupernodes; ++k )
        L.supernodes[k].workspace.Empty();
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
    for( int k=0; k<numSupernodes; ++k )
    {
        const symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[k];
        XSub.View( X, symbSN.offset, 0, symbSN.size, width );

        Matrix<F> d;
        L.supernodes[k].front.GetDiagonal( d );
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
  const numeric::LocalSymmFact<F>& L,
  F alpha, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLBackwardSolve");
#endif
    const int numSupernodes = S.supernodes.size();
    const int width = X.Width();
    for( int k=numSupernodes-1; k>=0; --k )
    {
        const symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[k];
        const numeric::LocalSymmFactSupernode<F>& numSN = L.supernodes[k];

        // Set up a workspace
        Matrix<F>& W = numSN.workspace;
        W.ResizeTo( numSN.front.Height(), width );
        Matrix<F> WT, WB;
        WT.View( W, 0, 0, symbSN.size, width );
        WB.View( W, symbSN.size, 0, W.Height()-symbSN.size, width );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.LockedView( X, symbSN.offset, 0, symbSN.size, X.Width() );
        WT = XT;

        // Update using the parent (if it exists)
        const int parentIndex = symbSN.parent;
        if( parentIndex != -1 )
        {
            Matrix<F>& parentWork = L.supernodes[parentIndex].workspace;
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
        }
        // else we are the root node

        // Call the custom supernode forward solve
        LocalSupernodeLDLBackwardSolve
        ( orientation, symbSN.size, alpha, numSN.front, W );

        // Store the supernode portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed
    for( int k=0; k<numSupernodes; ++k )
        L.supernodes[k].workspace.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<float>& L,
  float alpha, Matrix<float>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<float>& L,
  Matrix<float>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<float>& L,
  float alpha, Matrix<float>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<double>& L,
  double alpha, Matrix<double>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<double>& L,
  Matrix<double>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<double>& L,
  double alpha, Matrix<double>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& L,
  std::complex<float> alpha, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& L,
  Matrix<std::complex<float> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& L,
  std::complex<float> alpha, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& L,
  std::complex<double> alpha, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& L,
  Matrix<std::complex<double> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& L,
  std::complex<double> alpha, Matrix<std::complex<double> >& X );
