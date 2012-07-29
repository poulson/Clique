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

namespace cliq {

struct LocalSepOrLeaf
{
    int parent; // -1 if local root
    int offset;
    std::vector<int> indices;
};

struct DistSeparator
{
    mpi::Comm comm;
    int offset;
    std::vector<int> indices;
};

struct DistSeparatorTree
{
    // Full local binary tree
    //
    // NOTE: This is an array of pointers, as we will not know during 
    //       construction how many will need to be created
    std::vector<LocalSepOrLeaf*> localSepsAndLeaves;

    // One path through top of binary tree 
    //
    // NOTE: does not include the single-process separator/leaf
    std::vector<DistSeparator> distSeps;

    ~DistSeparatorTree()
    {
        if( std::uncaught_exception() )
        {
            std::cerr << "Uncaught exception in ~DistSepTree" << std::endl;
#ifndef RELEASE            
            DumpCallStack();
#endif
            return;
        }

        const int numLocal = localSepsAndLeaves.size();
        for( int i=0; i<numLocal; ++i )
            delete localSepsAndLeaves[i];

        const int numDist = distSeps.size();
        for( int i=0; i<numDist; ++i )
            mpi::CommFree( distSeps[i].comm );
    }
};

} // namespace cliq
