/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

struct SymmNodeInfo
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
    bool onLeft;
    int myOffset;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelIndices;
    // (maps from the child update indices to our frontal indices).
    std::vector<int> leftRelIndices, rightRelIndices;
};

struct FactorMetadata
{
    std::vector<int> numChildSendIndices;
    std::deque<int> leftColIndices, leftRowIndices,
                    rightColIndices, rightRowIndices;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::deque<int> > childRecvIndices;

    void EmptyChildRecvIndices() const
    {
        std::vector<std::deque<int> >().swap( childRecvIndices );
    }

    void Empty()
    {
        std::deque<int>().swap( leftColIndices );
        std::deque<int>().swap( leftRowIndices );
        std::deque<int>().swap( rightColIndices );
        std::deque<int>().swap( rightRowIndices );
        std::vector<int>().swap( numChildSendIndices );
        EmptyChildRecvIndices();
    }
};

struct SolveMetadata1d
{
    int localSize, localOffset;
    std::deque<int> leftIndices, rightIndices;
    std::vector<int> numChildSendIndices;
    std::vector<std::deque<int> > childRecvIndices;

    void Empty()
    {
        std::deque<int>().swap( leftIndices );
        std::deque<int>().swap( rightIndices );
        std::vector<int>().swap( numChildSendIndices );
        std::vector<std::deque<int> >().swap( childRecvIndices );
    }
};

struct SolveMetadata2d
{
    int localHeight, localWidth, localHeightOffset, localWidthOffset;
    std::deque<int> leftIndices, rightIndices;
    std::vector<int> numChildSendIndices;
    std::vector<std::deque<int> > childRecvIndices;

    void Empty()
    {
        std::deque<int>().swap( leftIndices );
        std::deque<int>().swap( rightIndices );
        std::vector<int>().swap( numChildSendIndices );
        std::vector<std::deque<int> >().swap( childRecvIndices );
    }
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
    int myOffset, leftSize, rightSize;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelIndices;

    // The relative indices of our child
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftRelIndices, rightRelIndices;

    FactorMetadata factorMeta;
    SolveMetadata1d solveMeta1d;
    SolveMetadata2d solveMeta2d;
};

struct DistSymmInfo
{
    std::vector<SymmNodeInfo> localNodes;
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
    CallStackEntry entry("DistSymmInfo::StoreReordered");
#endif
    const int localSize = distNodes.back().solveMeta1d.localOffset +
                          distNodes.back().solveMeta1d.localSize;
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
}

} // namespace cliq
