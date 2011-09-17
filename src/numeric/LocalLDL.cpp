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
void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalFactStruct& S, // can't be const due to map...
  numeric::LocalFactMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = S.lowerStructs.size();

    // Perform the local factorization
    for( int k=0; k<numSupernodes; ++k )
    {
        Matrix<F>& front = L.fronts[k];
        const int supernodeSize = S.sizes[k];
#ifndef RELEASE
        const int lowerStructSize = S.lowerStructs[k].size();
        if( front.Height() != supernodeSize+lowerStructSize ||
            front.Width()  != supernodeSize+lowerStructSize )
            throw std::logic_error("Front was not the proper size");
#endif

        // Add updates from children (if they exist)
        const int numChildren = S.children[k].size();
        if( numChildren == 2 )
        {
            const int leftIndex = S.children[k][0];
            const int rightIndex = S.children[k][1];
            const Matrix<F>& leftFront = L.fronts[leftIndex];
            const Matrix<F>& rightFront = L.fronts[rightIndex];
            const int leftSupernodeSize = S.sizes[leftIndex];
            const int rightSupernodeSize = S.sizes[rightIndex];
            const int leftUpdateSize = leftFront.Height()-leftSupernodeSize;
            const int rightUpdateSize = rightFront.Height()-rightSupernodeSize;
            
            const std::vector<int>& leftRelIndices = S.leftChildRelIndices[k];
            const std::vector<int>& rightRelIndices = S.rightChildRelIndices[k];

            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftFront, leftSupernodeSize, leftSupernodeSize,
              leftUpdateSize, leftUpdateSize );

            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightFront, rightSupernodeSize, rightSupernodeSize,
              rightUpdateSize, rightUpdateSize );

            // Add the left child's update matrix
            for( int jChild=0; jChild<leftUpdateSize; ++jChild )
            {
                const int jFront = leftRelIndices[jChild];
                for( int iChild=0; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = leftRelIndices[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
                    front.Update( iFront, jFront, -value );
                }
            }

            // Add the right child's update matrix
            for( int jChild=0; jChild<rightUpdateSize; ++jChild )
            {
                const int jFront = rightRelIndices[jChild];
                for( int iChild=0; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = rightRelIndices[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
                    front.Update( iFront, jFront, -value );
                }
            }
        }

        // Call the custom partial LDL
        LocalSupernodeLDL( orientation, front, supernodeSize );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalFactStruct& S,
  numeric::LocalFactMatrix<float>& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalFactStruct& S,
  numeric::LocalFactMatrix<double>& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalFactStruct& S,
  numeric::LocalFactMatrix<std::complex<float> >& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalFactStruct& S,
  numeric::LocalFactMatrix<std::complex<double> >& L );

