/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

// 'Supernode' should perhaps be preferred to 'node', but since we will always
// use supernodes, the extra verbage is unnecessarily cumbersome.

struct SymmNode
{
    int size, offset; 
    int parent; // -1 if root separator
    std::vector<int> children;
    std::vector<int> lowerStruct;
};

struct DistSymmNode
{
    bool onLeft; // irrelevant if root node
    mpi::Comm comm;
    int size, offset;
    std::vector<int> lowerStruct;
};

struct DistSymmElimTree
{
    // NOTE: This is an array of pointers, as we will not know how many 
    //       are needed during construction
    std::vector<SymmNode*> localNodes;
    std::vector<DistSymmNode> distNodes;

    ~DistSymmElimTree()
    {
        if( std::uncaught_exception() )
        {
            std::cerr << "Uncaught exception in ~DistSymmElimTree" << std::endl;
#ifndef RELEASE            
            DumpCallStack();
#endif
            return;
        }

        const int numLocal = localNodes.size();
        for( int i=0; i<numLocal; ++i )
            delete localNodes[i];

        const int numDist = distNodes.size();
        for( int i=0; i<numDist; ++i )
            mpi::CommFree( distNodes[i].comm );
    }
};

} // namespace cliq
