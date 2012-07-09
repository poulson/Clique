/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
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

namespace cliq {

void LocalSymmetricAnalysis( const DistSymmElimTree& eTree, DistSymmInfo& info )
{
#ifndef RELEASE
    PushCallStack("LocalSymmetricAnalysis");
#endif
    const int numNodes = eTree.localNodes.size();
    info.localNodes.resize( numNodes );

    // Perform the symbolic factorization
    int myOffset = 0;
    std::vector<int>::iterator it;
    std::vector<int> childrenStruct, partialStruct, fullStruct, nodeIndices;
    for( int s=0; s<numNodes; ++s )
    {
        const LocalSymmNode& node = eTree.localNodes[s];
        LocalSymmNodeInfo& nodeInfo = info.localNodes[s];
        nodeInfo.size = node.size;
        nodeInfo.offset = node.offset;
        nodeInfo.myOffset = myOffset;
        nodeInfo.parent = node.parent;
        nodeInfo.children = node.children;
        nodeInfo.origLowerStruct = node.lowerStruct;

        const int numChildren = node.children.size();
        const int numOrigLowerIndices = node.lowerStruct.size();
#ifndef RELEASE
        if( numChildren != 0 && numChildren != 2 )
            throw std::logic_error("Tree must be built from bisections");
#endif
        if( numChildren == 2 )
        {
            const int left = node.children[0];
            const int right = node.children[1];
            LocalSymmNodeInfo& leftChild = info.localNodes[left];
            LocalSymmNodeInfo& rightChild = info.localNodes[right];
            leftChild.isLeftChild = true;
            rightChild.isLeftChild = false;

            // Union the child lower structs
            const int numLeftIndices = leftChild.lowerStruct.size();
            const int numRightIndices = rightChild.lowerStruct.size();
#ifndef RELEASE
            for( int i=1; i<numLeftIndices; ++i )
            {
                if( leftChild.lowerStruct[i] <= leftChild.lowerStruct[i-1] )
                {
                    std::ostringstream msg;
                    msg << "Left child struct was not sorted for s=" << s 
                        << "\n";
                    throw std::logic_error( msg.str().c_str() );
                }
            }
            for( int i=1; i<numRightIndices; ++i )
            {
                if( rightChild.lowerStruct[i] <= rightChild.lowerStruct[i-1] )
                {
                    std::ostringstream msg;
                    msg << "Right child struct was not sorted for s=" << s
                        << "\n";
                    throw std::logic_error( msg.str().c_str() );
                }
            }
#endif
            const int childrenStructMaxSize = numLeftIndices + numRightIndices;
            childrenStruct.resize( childrenStructMaxSize );
            it = std::set_union
            ( leftChild.lowerStruct.begin(), leftChild.lowerStruct.end(),
              rightChild.lowerStruct.begin(), rightChild.lowerStruct.end(),
              childrenStruct.begin() );
            const int childrenStructSize = int(it-childrenStruct.begin());
            childrenStruct.resize( childrenStructSize );

            // Union the lower structure of this node
#ifndef RELEASE
            for( int i=1; i<numOrigLowerIndices; ++i )
            {
                if( node.lowerStruct[i] <= node.lowerStruct[i-1] )
                {
                    std::ostringstream msg;
                    msg << "Original struct not sorted for s=" << s << "\n";
                    throw std::logic_error( msg.str().c_str() );
                }
            }
#endif
            const int partialStructMaxSize = 
                childrenStructSize + numOrigLowerIndices;
            partialStruct.resize( partialStructMaxSize );
            it = std::set_union
            ( node.lowerStruct.begin(), node.lowerStruct.end(),
              childrenStruct.begin(), childrenStruct.end(),
              partialStruct.begin() );
            const int partialStructSize = int(it-partialStruct.begin());
            partialStruct.resize( partialStructSize );

            // Union again with the node indices
            nodeIndices.resize( node.size );
            for( int i=0; i<node.size; ++i )
                nodeIndices[i] = node.offset + i;
            fullStruct.resize( node.size + partialStructSize );
            it = std::set_union
            ( partialStruct.begin(), partialStruct.end(),
              nodeIndices.begin(), nodeIndices.end(),
              fullStruct.begin() );
            fullStruct.resize( int(it-fullStruct.begin()) );

            // Construct the relative indices of the original lower structure
            it = fullStruct.begin();
            nodeInfo.origLowerRelIndices.resize( numOrigLowerIndices );
            for( int i=0; i<numOrigLowerIndices; ++i )
            {
                const int index = node.lowerStruct[i];
                it = std::lower_bound ( it, fullStruct.end(), index );
                nodeInfo.origLowerRelIndices[i] = int(it-fullStruct.begin());
            }

            // Construct the relative indices of the children
            nodeInfo.leftChildRelIndices.resize( numLeftIndices );
            it = fullStruct.begin();
            for( int i=0; i<numLeftIndices; ++i )
            {
                const int index = leftChild.lowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                nodeInfo.leftChildRelIndices[i] = int(it-fullStruct.begin());
            }
            nodeInfo.rightChildRelIndices.resize( numRightIndices );
            it = fullStruct.begin();
            for( int i=0; i<numRightIndices; ++i )
            {
                const int index = rightChild.lowerStruct[i];
                it = std::lower_bound( it, fullStruct.end(), index );
                nodeInfo.rightChildRelIndices[i] = int(it-fullStruct.begin());
            }

            // Form lower struct of this node by removing node indices
            // (which take up the first node.size indices of fullStruct)
            const int lowerStructSize = fullStruct.size()-node.size;
            nodeInfo.lowerStruct.resize( lowerStructSize );
            for( int i=0; i<lowerStructSize; ++i )
                nodeInfo.lowerStruct[i] = fullStruct[node.size+i];
        }
        else // numChildren == 0, so this is a leaf node 
        {
            nodeInfo.lowerStruct = node.lowerStruct;
            
            // Construct the trivial relative indices of the original structure
            nodeInfo.origLowerRelIndices.resize( numOrigLowerIndices );
            for( int i=0; i<numOrigLowerIndices; ++i )
                nodeInfo.origLowerRelIndices[i] = i + nodeInfo.size;
        }

        myOffset += nodeInfo.size;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
