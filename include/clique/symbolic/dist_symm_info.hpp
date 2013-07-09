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
    std::vector<int> origLowerRelInd;
    // (maps from the child update indices to our frontal indices).
    std::vector<int> leftRelInd, rightRelInd;
};

struct FactorMetadata
{
    std::vector<int> numChildSendInd;
    std::deque<int> leftColInd, leftRowInd, rightColInd, rightRowInd;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::deque<int> > childRecvInd;

    void EmptyChildRecvIndices() const
    {
        std::vector<std::deque<int> >().swap( childRecvInd );
    }

    void Empty()
    {
        std::deque<int>().swap( leftColInd );
        std::deque<int>().swap( leftRowInd );
        std::deque<int>().swap( rightColInd );
        std::deque<int>().swap( rightRowInd );
        std::vector<int>().swap( numChildSendInd );
        EmptyChildRecvIndices();
    }
};

struct SolveMetadata1d
{
    int localSize;
    std::deque<int> leftInd, rightInd;
    std::vector<int> numChildSendInd;
    std::vector<std::deque<int> > childRecvInd;

    void Empty()
    {
        std::deque<int>().swap( leftInd );
        std::deque<int>().swap( rightInd );
        std::vector<int>().swap( numChildSendInd );
        std::vector<std::deque<int> >().swap( childRecvInd );
    }
};

struct SolveMetadata2d
{
    int localHeight, localWidth;
    std::deque<int> leftInd, rightInd;
    std::vector<int> numChildSendInd;
    std::vector<std::deque<int> > childRecvInd;

    void Empty()
    {
        std::deque<int>().swap( leftInd );
        std::deque<int>().swap( rightInd );
        std::vector<int>().swap( numChildSendInd );
        std::vector<std::deque<int> >().swap( childRecvInd );
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
    std::vector<int> origLowerRelInd;

    // The relative indices of our child
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftRelInd, rightRelInd;

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
    int localSize=0;
    const int numLocal = localNodes.size();
    for( int s=0; s<numLocal; ++s )
        localSize += localNodes[s].size;
    const int numDist = distNodes.size();
    for( int s=1; s<numDist; ++s )
        localSize += distNodes[s].solveMeta1d.localSize;
    reordered.resize( localSize );

    int localOffset=0;
    for( int s=0; s<numLocal; ++s )
    {
        const int size = localNodes[s].size;
        const int offset = localNodes[s].offset;
        for( int j=0; j<size; ++j )
            reordered[localOffset++] = j+offset;
    }
    for( int s=1; s<numDist; ++s )
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
