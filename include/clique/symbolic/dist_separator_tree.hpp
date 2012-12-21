/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
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
