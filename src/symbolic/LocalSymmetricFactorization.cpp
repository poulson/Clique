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

void clique::symbolic::LocalSymmetricFactorization
( const LocalOrigStruct& localOrig,
        LocalFactStruct& localFact )
{
#ifndef RELEASE
    PushCallStack("symbolic::LocalSymmetricFactorization");
#endif
    const int numSupernodes = localOrig.sizes.size();
    localFact.sizes = localOrig.sizes;
    localFact.offsets = localOrig.sizes;
    localFact.children = localOrig.children;
    localFact.parents = localOrig.parents;

    localFact.lowerStructs.resize( numSupernodes );
    localFact.origLowerRelIndices.resize( numSupernodes );
    localFact.leftChildRelIndices.resize( numSupernodes );
    localFact.rightChildRelIndices.resize( numSupernodes );
    localFact.isLeftChild.resize( numSupernodes );

    // Perform the local symbolic factorization
    std::vector<int>::iterator it;
    std::vector<int> childrenStruct, partialStruct, fullStruct,
                     supernodeIndices;
    for( int k=0; k<numSupernodes; ++k )
    {
        const int supernodeSize = localOrig.sizes[k];
        const int supernodeOffset = localOrig.offsets[k];
        const int numChildren = localOrig.children[k].size();
        const std::vector<int>& origLowerStruct = localOrig.lowerStructs[k];
        std::vector<int>& lowerStruct = localFact.lowerStructs[k];

#ifndef RELEASE
        if( numChildren != 0 && numChildren != 2 )
            throw std::logic_error("Tree must be built from bisections");
#endif

        if( numChildren == 2 )
        {
            const int leftIndex = localOrig.children[k][0];
            const int rightIndex = localOrig.children[k][1];
            localFact.isLeftChild[leftIndex] = true;
            localFact.isLeftChild[rightIndex] = false;

            const std::vector<int>& leftChildLowerStruct = 
                localFact.lowerStructs[leftIndex];
            const std::vector<int>& rightChildLowerStruct = 
                localFact.lowerStructs[rightIndex];
            const int numLeftIndices = leftChildLowerStruct.size();
            const int numRightIndices = rightChildLowerStruct.size();

            // Union the child lower structs
            const int childrenStructMaxSize = numLeftIndices + numRightIndices;
            childrenStruct.resize( childrenStructMaxSize );
            it = std::set_union
            ( leftChildLowerStruct.begin(), leftChildLowerStruct.end(),
              rightChildLowerStruct.begin(), rightChildLowerStruct.end(),
              childrenStruct.begin() );
            childrenStruct.resize( int(it-childrenStruct.begin()) );

            // Union the lower structure of this supernode
            const int partialStructMaxSize = 
                childrenStruct.size() + origLowerStruct.size();
            partialStruct.resize( partialStructMaxSize );
            it = std::set_union
            ( origLowerStruct.begin(), origLowerStruct.end(),
              childrenStruct.begin(), childrenStruct.end(),
              partialStruct.begin() );
            partialStruct.resize( int(it-partialStruct.begin()) );

            // Union again with the supernode indices
            supernodeIndices.resize( supernodeSize );
            for( int i=0; i<supernodeSize; ++i )
                supernodeIndices[i] = supernodeOffset + i;
            fullStruct.resize( supernodeSize + partialStruct.size() );
            it = std::set_union
            ( partialStruct.begin(), partialStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              fullStruct.begin() );
            fullStruct.resize( int(it-fullStruct.begin()) );

            // Construct the relative indices of the original lower structure
            const int numOrigLowerIndices = origLowerStruct.size();
            it = fullStruct.begin();
            for( int i=0; i<numOrigLowerIndices; ++i )
            {
                const int index = origLowerStruct[i];
                it = std::lower_bound ( it, fullStruct.end(), index );
                localFact.origLowerRelIndices[k][index] = *it;
            }

            // Construct the relative indices of the children
            localFact.leftChildRelIndices[k].resize( numLeftIndices );
            it = fullStruct.begin();
            for( int i=0; i<numLeftIndices; ++i )
            {
                const int index = leftChildLowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                localFact.leftChildRelIndices[k][i] = *it;
            }
            localFact.rightChildRelIndices[k].resize( numRightIndices );
            it = fullStruct.begin();
            for( int i=0; i<numRightIndices; ++it )
            {
                const int index = rightChildLowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                localFact.rightChildRelIndices[k][i] = *it;
            }

            // Form lower struct of this supernode by removing supernode indices
            lowerStruct.resize( fullStruct.size() );
            it = std::set_difference
            ( fullStruct.begin(), fullStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              lowerStruct.begin() );
            lowerStruct.resize( int(it-lowerStruct.begin()) );
        }
        else // numChildren == 0, so this is a leaf supernode 
        {
            localFact.lowerStructs[k] = localOrig.lowerStructs[k];        
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

