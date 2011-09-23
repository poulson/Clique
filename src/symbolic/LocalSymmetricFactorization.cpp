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
( const LocalSymmOrig& orig, LocalSymmFact& fact )
{
#ifndef RELEASE
    PushCallStack("symbolic::LocalSymmetricFactorization");
#endif
    const int numSupernodes = orig.supernodes.size();
    fact.supernodes.resize( numSupernodes );

    // Perform the symbolic factorization
    int myOffset = 0;
    std::vector<int>::iterator it;
    std::vector<int> childrenStruct, partialStruct, fullStruct,
                     supernodeIndices;
    for( int k=0; k<numSupernodes; ++k )
    {
        const LocalSymmOrigSupernode& origSN = orig.supernodes[k];
        LocalSymmFactSupernode& factSN = fact.supernodes[k];
        factSN.size = origSN.size;
        factSN.offset = origSN.offset;
        factSN.myOffset = myOffset;
        factSN.parent = origSN.parent;
        factSN.children = origSN.children;

        const int numChildren = origSN.children.size();
#ifndef RELEASE
        if( numChildren != 0 && numChildren != 2 )
            throw std::logic_error("Tree must be built from bisections");
#endif
        if( numChildren == 2 )
        {
            const int leftIndex = origSN.children[0];
            const int rightIndex = origSN.children[1];
            LocalSymmFactSupernode& leftChild = fact.supernodes[leftIndex];
            LocalSymmFactSupernode& rightChild = fact.supernodes[rightIndex];
            leftChild.isLeftChild = true;
            rightChild.isLeftChild = false;

            // Union the child lower structs
            const int numLeftIndices = leftChild.lowerStruct.size();
            const int numRightIndices = rightChild.lowerStruct.size();
            const int childrenStructMaxSize = numLeftIndices + numRightIndices;
            childrenStruct.resize( childrenStructMaxSize );
            it = std::set_union
            ( leftChild.lowerStruct.begin(), leftChild.lowerStruct.end(),
              rightChild.lowerStruct.begin(), rightChild.lowerStruct.end(),
              childrenStruct.begin() );
            const int childrenStructSize = int(it-childrenStruct.begin());
            childrenStruct.resize( childrenStructSize );

            // Union the lower structure of this supernode
            const int numOrigLowerIndices = origSN.lowerStruct.size();
            const int partialStructMaxSize = 
                childrenStructSize + numOrigLowerIndices;
            partialStruct.resize( partialStructMaxSize );
            it = std::set_union
            ( origSN.lowerStruct.begin(), origSN.lowerStruct.end(),
              childrenStruct.begin(), childrenStruct.end(),
              partialStruct.begin() );
            const int partialStructSize = int(it-partialStruct.begin());
            partialStruct.resize( partialStructSize );

            // Union again with the supernode indices
            supernodeIndices.resize( origSN.size );
            for( int i=0; i<origSN.size; ++i )
                supernodeIndices[i] = origSN.offset + i;
            fullStruct.resize( origSN.size + partialStructSize );
            it = std::set_union
            ( partialStruct.begin(), partialStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              fullStruct.begin() );
            fullStruct.resize( int(it-fullStruct.begin()) );

            // Construct the relative indices of the original lower structure
            it = fullStruct.begin();
            for( int i=0; i<numOrigLowerIndices; ++i )
            {
                const int index = origSN.lowerStruct[i];
                it = std::lower_bound ( it, fullStruct.end(), index );
                factSN.origLowerRelIndices[index] = *it;
            }

            // Construct the relative indices of the children
            factSN.leftChildRelIndices.resize( numLeftIndices );
            it = fullStruct.begin();
            for( int i=0; i<numLeftIndices; ++i )
            {
                const int index = leftChild.lowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                factSN.leftChildRelIndices[i] = *it;
            }
            factSN.rightChildRelIndices.resize( numRightIndices );
            it = fullStruct.begin();
            for( int i=0; i<numRightIndices; ++i )
            {
                const int index = rightChild.lowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                factSN.rightChildRelIndices[i] = *it;
            }

            // Form lower struct of this supernode by removing supernode indices
            factSN.lowerStruct.resize( fullStruct.size() );
            it = std::set_difference
            ( fullStruct.begin(), fullStruct.end(),
              supernodeIndices.begin(), supernodeIndices.end(),
              factSN.lowerStruct.begin() );
            factSN.lowerStruct.resize( int(it-factSN.lowerStruct.begin()) );
        }
        else // numChildren == 0, so this is a leaf supernode 
        {
            factSN.lowerStruct = origSN.lowerStruct;
        }

        myOffset += factSN.size;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

