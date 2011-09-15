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
        symbolic::LocalFactStruct& SLocal, // can't be const due to map...
  const numeric::LocalOrigMatrix<F>& ALocal,
        numeric::LocalFactMatrix<F>& LLocal )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = SLocal.lowerStructs.size();

    // Perform the local factorization
    for( int k=0; k<numSupernodes; ++k )
    {
        const int supernodeOffset = SLocal.offsets[k];
        const int supernodeSize = SLocal.sizes[k];
        const int lowerStructSize = SLocal.lowerStructs[k].size();

        const std::vector<F>& nonzeros = ALocal.nonzeros[k];
        const std::vector<int>& colOffsets = ALocal.colOffsets[k];
        const std::vector<int>& rowIndices = ALocal.rowIndices[k];
        std::map<int,int>& origLowerRelIndices = SLocal.origLowerRelIndices[k];
        const std::vector<int>& leftChildRelIndices = 
            SLocal.leftChildRelIndices[k];
        const std::vector<int>& rightChildRelIndices = 
            SLocal.rightChildRelIndices[k];

        // Expand the original sparse matrix into the frontal matrix.
        Matrix<F>& front = LLocal.fronts[k];
        front.ResizeTo
        ( supernodeSize+lowerStructSize, supernodeSize+lowerStructSize );
        front.SetToZero();
        const int numColumns = colOffsets.size()-1;
        for( int jFront=0; jFront<numColumns; ++jFront )
        {
            const int thisColOffset = colOffsets[jFront];
            const int thisColSize = colOffsets[jFront+1]-thisColOffset;
            const int* theseRowIndices = &rowIndices[thisColOffset];
            const F* theseNonzeros = &nonzeros[k];

            for( int iPacked=0; iPacked<thisColSize; ++iPacked )
            {
                const int iGlobal = theseRowIndices[iPacked];
                const int iFront = 
                    ( iGlobal < supernodeOffset+supernodeSize ? 
                      iGlobal : origLowerRelIndices[iGlobal] );

                front.Set( iFront, jFront, theseNonzeros[iPacked] );
            }
        }

        // Add updates from children (if they exist)
        const int numChildren = SLocal.children[k].size();
        if( numChildren == 2 )
        {
            const int leftChildIndex = SLocal.children[k][0];
            const int rightChildIndex = SLocal.children[k][1];
            const Matrix<F>& leftChildFront = LLocal.fronts[leftChildIndex];
            const Matrix<F>& rightChildFront = LLocal.fronts[rightChildIndex];
            const int leftChildSupernodeSize = SLocal.sizes[leftChildIndex];
            const int rightChildSupernodeSize = SLocal.sizes[rightChildIndex];
            const int leftChildUpdateSize = 
                leftChildFront.Height()-leftChildSupernodeSize;
            const int rightChildUpdateSize = 
                rightChildFront.Height()-rightChildSupernodeSize;

            Matrix<F> leftChildUpdate;
            leftChildUpdate.LockedView
            ( leftChildFront, leftChildSupernodeSize, leftChildSupernodeSize,
              leftChildUpdateSize, leftChildUpdateSize );

            Matrix<F> rightChildUpdate;
            rightChildUpdate.LockedView
            ( rightChildFront, rightChildSupernodeSize, rightChildSupernodeSize,
              rightChildUpdateSize, rightChildUpdateSize );

            // Add the left child's update matrix
            for( int jChild=0; jChild<leftChildUpdateSize; ++jChild )
            {
                const int jFront = leftChildRelIndices[jChild];
                for( int iChild=0; iChild<leftChildUpdateSize; ++iChild )
                {
                    const int iFront = leftChildRelIndices[iChild];
                    const F value = leftChildUpdate.Get(iChild,jChild);
                    front.Update( iFront, jFront, -value );
                }
            }

            // Add the right child's update matrix
            for( int jChild=0; jChild<rightChildUpdateSize; ++jChild )
            {
                const int jFront = rightChildRelIndices[jChild];
                for( int iChild=0; iChild<rightChildUpdateSize; ++iChild )
                {
                    const int iFront = rightChildRelIndices[iChild];
                    const F value = rightChildUpdate.Get(iChild,jChild);
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
        symbolic::LocalFactStruct& SLocal,
  const numeric::LocalOrigMatrix<float>& ALocal,
        numeric::LocalFactMatrix<float>& LLocal );

template void clique::numeric::LocalLDL
( Orientation orientation,
        symbolic::LocalFactStruct& SLocal,
  const numeric::LocalOrigMatrix<double>& ALocal,
        numeric::LocalFactMatrix<double>& LLocal );

template void clique::numeric::LocalLDL
( Orientation orientation,
        symbolic::LocalFactStruct& SLocal,
  const numeric::LocalOrigMatrix<std::complex<float> >& ALocal,
        numeric::LocalFactMatrix<std::complex<float> >& LLocal );

template void clique::numeric::LocalLDL
( Orientation orientation,
        symbolic::LocalFactStruct& SLocal,
  const numeric::LocalOrigMatrix<std::complex<double> >& ALocal,
        numeric::LocalFactMatrix<std::complex<double> >& LLocal );

