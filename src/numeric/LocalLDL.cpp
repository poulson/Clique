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
  const symbolic::LocalFactStruct& SLocal,
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
        // Expand the original sparse matrix into the frontal matrix.
        const int supernodeOffset = SLocal.offsets[k];
        const int supernodeSize = SLocal.sizes[k];
        const int lowerStructSize = SLocal.lowerStructs[k].size();
        LLocal.fronts[k].ResizeTo( supernodeSize + lowerStructSize );
        LLocal.fronts[k].SetToZero();
        const int numColumns = ALocal.colOffsets[k].size()-1;
        for( int jPacked=0; jPacked<numColumns; ++jPacked )
        {
            const int colOffset = ALocal.colOffsets[k][jPacked];
            const int numRowIndices = ALocal.colOffsets[k][jPacked+1]-colOffset;
            const int* rowIndices = &ALocal.rowIndices[k][colOffset];

            // HERE

            for( int iPacked=0; iPacked<numRowIndices; ++iPacked )
            {
                // HERE
            }
        }

        // Add updates from children (if they exist)
        const int numChildren = SLocal.children[k].size();

        // Call the custom partial LDL
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

