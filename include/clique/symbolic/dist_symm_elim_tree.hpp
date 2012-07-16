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
#ifndef CLIQUE_DIST_SYMM_ELIM_TREE_HPP
#define CLIQUE_DIST_SYMM_ELIM_TREE_HPP 1

namespace cliq {

// 'Supernode' should perhaps be preferred to 'node', but since we will always
// use supernodes, the extra verbage is unnecessarily cumbersome.

struct LocalSymmNode
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
    std::vector<LocalSymmNode*> localNodes;
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

#endif /* CLIQUE_DIST_SYMM_ELIM_TREE_HPP */
