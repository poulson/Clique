/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
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

void clique::symbolic::BottomSymmetricFactorization
( const BottomOrigStruct& bottomOrig,
        BottomFactStruct& bottomFact )
{
#ifndef RELEASE
    PushCallStack("symbolic::BottomSymmetricFactorization");
#endif
    const int numSupernodes = bottomOrig.sizes.size();
    bottomFact.sizes = bottomOrig.sizes;
    bottomFact.offsets = bottomOrig.sizes;
    bottomFact.lowerStructs.resize( numSupernodes );
    bottomFact.children.resize( numSupernodes );
    bottomFact.leftChildRelIndices.resize( numSupernodes );
    bottomFact.rightChildRelIndices.resize( numSupernodes );

    // Perform the local symbolic factorization
    std::vector<int> childrenStruct, partialStruct, fullStruct,
                     supernodeIndices;
    for( int k=0; k<numSupernodes; ++k )
    {
        const int supernodeSize = bottomOrig.sizes[k];
        const int supernodeOffset = bottomOrig.offsets[k];
        const int numChildren = bottomOrig.children[k].size();
        const std::vector<int>& origLowerStruct = bottomOrig.lowerStructs[k];
        std::vector<int>& lowerStruct = bottomFact.lowerStructs[k];

        bottomFact.children[k] = bottomOrig.children[k];
#ifndef RELEASE
        if( numChildren != 0 && numChildren != 2 )
            throw std::logic_error("Tree must be built from bisections");
#endif

        if( numChildren == 2 )
        {
            const std::vector<int>& leftChildLowerStruct = 
                bottomFact.lowerStructs[bottomOrig.children[k][0]];
            const std::vector<int>& rightChildLowerStruct = 
                bottomFact.lowerStructs[bottomOrig.children[k][1]];
            const int numLeftIndices = leftChildLowerStruct.size();
            const int numRightIndices = rightChildLowerStruct.size();

            // Union the child lower structs
            const int childrenStructMaxSize = numLeftIndices + numRightIndices;
            childrenStruct.resize( childrenStructMaxSize );
            std::set_union
            ( leftChildLowerStruct.begin(), leftChildLowerStruct.end(),
              rightChildLowerStruct.begin(), rightChildLowerStruct.end(),
              childrenStruct.begin() );

            // Union the lower structure of this supernode
            const int partialStructMaxSize = 
                childrenStruct.size() + origLowerStruct.size();
            partialStruct.resize( partialStructMaxSize );
            std::set_union
            ( origLowerStruct.begin(), origLowerStruct.end(),
              childrenStruct.begin(), childrenStruct.end(),
              partialStruct.begin() );

            // Union again with the supernode indices
            supernodeIndices.resize( supernodeSize );
            for( int i=0; i<supernodeSize; ++i )
                supernodeIndices[i] = supernodeOffset + i;
            fullStruct.resize( supernodeSize + partialStruct.size() );
            std::set_union
            ( partialStruct.begin(), partialStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              fullStruct.begin() );

            // Construct the relative indices of the children
            bottomFact.leftChildRelIndices[k].resize( numLeftIndices );
            bottomFact.rightChildRelIndices[k].resize( numRightIndices );
            std::vector<int>::iterator it = fullStruct.begin();
            for( int i=0; i<numLeftIndices; ++i )
            {
                it = std::lower_bound
                     ( it, fullStruct.end(), leftChildLowerStruct[i] );
                bottomFact.leftChildRelIndices[k][i] = *it;
            }
            it = fullStruct.begin();
            for( int i=0; i<numRightIndices; ++it )
            {
                it = std::lower_bound
                     ( it, fullStruct.end(), rightChildLowerStruct[i] );
                bottomFact.rightChildRelIndices[k][i] = *it;
            }

            // Form lower struct of this supernode by removing supernode indices
            lowerStruct.resize( fullStruct.size() );
            std::set_difference
            ( fullStruct.begin(), fullStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              lowerStruct.begin() );
        }
        else // numChildren == 0, so this is a leaf supernode 
        {
            bottomFact.lowerStructs[k] = bottomOrig.lowerStructs[k];        
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

