/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

struct LocalSymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, offset; 
    int parent; // -1 if root separator
    std::vector<int> children;
    std::vector<int> origLowerStruct;

    //
    // The following is computed during analysis
    //
    bool isLeftChild;
    int myOffset;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelIndices;
    // (maps from the child update indices to our frontal indices).
    std::vector<int> leftChildRelIndices, rightChildRelIndices;
};

struct DistSymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, offset;
    std::vector<int> origLowerStruct;
    bool onLeft;
    mpi::Comm comm;

    //
    // The following is computed during analysis
    //
    Grid* grid;
    int myOffset, leftChildSize, rightChildSize;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelIndices;

    // The relative indices of our child
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftChildRelIndices, rightChildRelIndices;

    // Helpers for the factorization
    std::vector<int> numChildFactSendIndices;
    std::deque<int> leftChildFactColIndices, leftChildFactRowIndices,
                    rightChildFactColIndices, rightChildFactRowIndices;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::deque<int> > childFactRecvIndices;

    // Helpers for solving with 1d right-hand sides
    std::deque<int> leftChildSolveIndices, rightChildSolveIndices;
    int localSize1d, localOffset1d;
    std::vector<int> numChildSolveSendIndices;
    std::vector<std::deque<int> > childSolveRecvIndices;
};

struct DistSymmInfo
{
    std::vector<LocalSymmNodeInfo> localNodes;
    std::vector<DistSymmNodeInfo> distNodes;
    ~DistSymmInfo();

    void StoreReordered( std::vector<int>& reordered ) const;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline
DistSymmInfo::~DistSymmInfo()
{
    const int numDist = distNodes.size();
    for( int s=0; s<numDist; ++s )
    {
        delete distNodes[s].grid;
        mpi::CommFree( distNodes[s].comm );
    }
}

inline void
DistSymmInfo::StoreReordered( std::vector<int>& reordered ) const
{
#ifndef RELEASE
    PushCallStack("DistSymmInfo::StoreReordered");
#endif
    const int localSize = distNodes.back().localOffset1d +
                          distNodes.back().localSize1d;
    reordered.resize( localSize );

    int localOffset=0;
    const int numDistNodes = distNodes.size();
    const int numLocalNodes = localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const int size = localNodes[s].size;
        const int offset = localNodes[s].offset;
        for( int j=0; j<size; ++j )
            reordered[localOffset++] = j+offset;
    }
    for( int s=1; s<numDistNodes; ++s )
    {
        const int size = distNodes[s].size;
        const int offset = distNodes[s].offset;
        const Grid& grid = *distNodes[s].grid;
        const int gridSize = grid.Size();
        const int gridRank = grid.Rank();
        // ASSUMPTION: The alignment is zero
        for( int j=gridRank; j<size; j+=gridSize )
            reordered[localOffset++] = j+offset;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
