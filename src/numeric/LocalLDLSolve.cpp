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
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<F>& L,
  F alpha, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
#endif
    const int numSupernodes = S.lowerStructs.size();
    const int width = X.Width();

    L.solutions.resize( numSupernodes );
    for( int k=0; k<numSupernodes; ++k )
    {
        const int numChildren = S.children[k].size();
        const int supernodeSize = S.sizes[k];
        const int supernodeOffset = S.offsets[k];
        const Matrix<F>& front = L.fronts[k];
        Matrix<F>& Y = L.solutions[k];

        Matrix<F> XT;
        XT.LockedView( X, supernodeOffset, 0, supernodeSize, X.Width() );
        
        Y.ResizeTo( front.Height(), width );
        Matrix<F> YT, YB;
        YT.View( Y, 0, 0, supernodeSize, width );
        YB.View( Y, supernodeSize, 0, Y.Height()-supernodeSize, width );

        // Pull in the relevant information from the RHS
        YT = XT;
        YB.SetToZero();

        // Update using the children (if they exist)
        if( numChildren == 2 )
        {
            const int leftIndex = S.children[k][0];
            const int rightIndex = S.children[k][1];
            Matrix<F>& leftSol = L.solutions[leftIndex];
            Matrix<F>& rightSol = L.solutions[rightIndex];
            const int leftSupernodeSize = S.sizes[leftIndex];
            const int rightSupernodeSize = S.sizes[rightIndex];
            const int leftUpdateSize = leftSol.Height()-leftSupernodeSize;
            const int rightUpdateSize = rightSol.Height()-rightSupernodeSize;

            const std::vector<int>& leftRelIndices = S.leftChildRelIndices[k];
            const std::vector<int>& rightRelIndices = S.rightChildRelIndices[k];

            // Add the left child's portion of the solution onto ours
            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftSol, leftSupernodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = leftRelIndices[iChild]; 
                for( int j=0; j<width; ++j )
                    Y.Update( iFront, j, -leftUpdate.Get(iFront,j) );
            }
            leftSol.Empty();

            // Add the right child's portion of the solution onto ours
            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightSol, rightSupernodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = rightRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    Y.Update( iFront, j, -rightUpdate.Get(iFront,j) );
            }
            rightSol.Empty();
        }
        // else numChildren == 0

        // Call the custom supernode forward solve
        LocalSupernodeLDLForwardSolve( supernodeSize, alpha, front, Y );

        // Store the supernode portion of the result
        XT = YT;
    }
    L.solutions.clear();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<F>& L,
  Matrix<F>& X, bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLDiagonalSolve");
#endif
    const int numSupernodes = S.lowerStructs.size();
    const int width = X.Width();

    Matrix<F> XSub;
    for( int k=0; k<numSupernodes; ++k )
    {
        const int supernodeSize = S.sizes[k];
        const int supernodeOffset = S.offsets[k];
        XSub.View( X, supernodeOffset, 0, supernodeSize, width );

        // Matrix<F> d;
        // L.fronts[k].GetDiagonal( d );
        // LocalSupernodeLDLDiagonalSolve
        // ( supernodeSize, d, XSub, checkIfSingular );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalFactStruct& S, 
  const numeric::LocalFactMatrix<F>& L,
  F alpha, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLBackwardSolve");
#endif
    const int numSupernodes = S.lowerStructs.size();
    const int width = X.Width();

    L.solutions.resize( numSupernodes );
    for( int k=numSupernodes-1; k>=0; --k )
    {
        const int numChildren = S.children[k].size();
        const int supernodeSize = S.sizes[k];
        const int supernodeOffset = S.offsets[k];
        const Matrix<F>& front = L.fronts[k];
        Matrix<F>& Y = L.solutions[k];

        Matrix<F> XT;
        XT.LockedView( X, supernodeOffset, 0, supernodeSize, X.Width() );
        
        Y.ResizeTo( front.Height(), width );
        Matrix<F> YT, YB;
        YT.View( Y, 0, 0, supernodeSize, width );
        YB.View( Y, supernodeSize, 0, Y.Height()-supernodeSize, width );

        // Pull in the relevant information from the RHS
        YT = XT;

        // Update using the parent (if it exists)
        const int parentIndex = S.parents[k];
        if( parentIndex != -1 )
        {
            Matrix<F>& parentSol = L.solutions[parentIndex];
            const bool isLeftChild = S.isLeftChild[k];
            const int parentSupernodeSize = S.sizes[parentIndex];
            const int parentUpdateSize = parentSol.Height()-parentSupernodeSize;
            const int currentUpdateSize = YB.Height();

            const std::vector<int>& parentRelIndices = 
              ( isLeftChild ? 
                S.leftChildRelIndices[parentIndex] :
                S.rightChildRelIndices[parentIndex] );

            for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
            {
                const int iParent = parentRelIndices[iCurrent];
                for( int j=0; j<width; ++j )
                    YB.Set( iCurrent, j, parentSol.Get(iParent,j) );
            }

            // The left child is numbered lower than the right child, so 
            // we can safely free the parent's solution if we are the left child
            if( isLeftChild )
                parentSol.Empty();
        }
        // else we are the root node

        // Call the custom supernode forward solve
        LocalSupernodeLDLBackwardSolve
        ( orientation, supernodeSize, alpha, front, Y );

        // Store the supernode portion of the result
        XT = YT;
    }
    L.solutions.clear();
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<float>& L,
  float alpha, Matrix<float>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<float>& L,
  Matrix<float>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<float>& L,
  float alpha, Matrix<float>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<double>& L,
  double alpha, Matrix<double>& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<double>& L,
  Matrix<double>& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<double>& L,
  double alpha, Matrix<double>& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<std::complex<float> >& L,
  std::complex<float> alpha, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<std::complex<float> >& L,
  Matrix<std::complex<float> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<std::complex<float> >& L,
  std::complex<float> alpha, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalLDLForwardSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<std::complex<double> >& L,
  std::complex<double> alpha, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalLDLDiagonalSolve
( const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<std::complex<double> >& L,
  Matrix<std::complex<double> >& X, bool checkIfSingular );
template void clique::numeric::LocalLDLBackwardSolve
( Orientation orientation,
  const symbolic::LocalFactStruct& S,
  const numeric::LocalFactMatrix<std::complex<double> >& L,
  std::complex<double> alpha, Matrix<std::complex<double> >& X );
