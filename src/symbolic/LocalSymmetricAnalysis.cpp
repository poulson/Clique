/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "clique.hpp"

namespace cliq {

void LocalSymmetricAnalysis( const DistSymmElimTree& eTree, DistSymmInfo& info )
{
#ifndef RELEASE
    CallStackEntry entry("LocalSymmetricAnalysis");
#endif
    const int numNodes = eTree.localNodes.size();
    info.localNodes.resize( numNodes );

    // Perform the symbolic factorization
    int myOffset = 0;
    for( int s=0; s<numNodes; ++s )
    {
        const SymmNode& node = *eTree.localNodes[s];
        SymmNodeInfo& nodeInfo = info.localNodes[s];
        nodeInfo.size = node.size;
        nodeInfo.offset = node.offset;
        nodeInfo.myOffset = myOffset;
        nodeInfo.parent = node.parent;
        nodeInfo.children = node.children;
        nodeInfo.origLowerStruct = node.lowerStruct;

        const int numChildren = node.children.size();
#ifndef RELEASE
        if( numChildren != 0 && numChildren != 2 )
            throw std::logic_error("Tree must be built from bisections");
#endif
        if( numChildren == 2 )
        {
            const int left = node.children[0];
            const int right = node.children[1];
            SymmNodeInfo& leftChild = info.localNodes[left];
            SymmNodeInfo& rightChild = info.localNodes[right];
            leftChild.onLeft = true;
            rightChild.onLeft = false;
#ifndef RELEASE
            if( !IsStrictlySorted(leftChild.lowerStruct) )
            {
                if( IsSorted(leftChild.lowerStruct) )
                    throw std::logic_error("Repeat in left lower struct");
                else
                    throw std::logic_error("Left lower struct not sorted");
            }
            if( !IsStrictlySorted(rightChild.lowerStruct) )
            {
                if( IsSorted(rightChild.lowerStruct) )
                    throw std::logic_error("Repeat in right lower struct");
                else
                    throw std::logic_error("Right lower struct not sorted");
            }
            if( !IsStrictlySorted(node.lowerStruct) )
            {
                if( IsSorted(node.lowerStruct) )
                    throw std::logic_error("Repeat in original lower struct");
                else
                    throw std::logic_error("Original lower struct not sorted");
            }
#endif

            // Combine the structures of the children
            std::vector<int> childrenStruct;
            Union
            ( childrenStruct, leftChild.lowerStruct, rightChild.lowerStruct );

            // Now add in the original lower structure
            std::vector<int> partialStruct;
            Union( partialStruct, node.lowerStruct, childrenStruct );

            // Now the node indices
            std::vector<int> nodeIndices( node.size );
            for( int i=0; i<node.size; ++i )
                nodeIndices[i] = node.offset + i;
            std::vector<int> fullStruct;
            Union( fullStruct, partialStruct, nodeIndices );

            // Construct the relative indices of the original lower structure
            RelativeIndices
            ( nodeInfo.origLowerRelIndices, node.lowerStruct, fullStruct );

            // Construct the relative indices of the children
            RelativeIndices
            ( nodeInfo.leftRelIndices, leftChild.lowerStruct, fullStruct );
            RelativeIndices
            ( nodeInfo.rightRelIndices, rightChild.lowerStruct, fullStruct );

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
            const int numOrigLowerIndices = node.lowerStruct.size();
            nodeInfo.origLowerRelIndices.resize( numOrigLowerIndices );
            for( int i=0; i<numOrigLowerIndices; ++i )
                nodeInfo.origLowerRelIndices[i] = i + nodeInfo.size;
        }

        myOffset += nodeInfo.size;
    }
}

} // namespace cliq
