/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix()
: height_(0), width_(0)
{ }

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::DistNodalMatrix");
#endif
    Pull( inverseMap, info, X );
}

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix( const DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::DistNodalMatrix");
#endif
    *this = X;
    commMetas.clear();
}

template<typename F>
inline const DistNodalMatrix<F>&
DistNodalMatrix<F>::operator=( const DistNodalMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::operator=");
#endif
    commMetas.clear();
    height_ = X.Height();
    width_ = X.Width();

    // Copy over the nontrivial distributed nodes
    const int numDist = X.distNodes.size();
    distNodes.resize( numDist );
    for( int s=0; s<numDist; ++s )
    {
        distNodes[s].SetGrid( X.distNodes[s].Grid() );
        distNodes[s] = X.distNodes[s];
    }

    // Copy over the local nodes
    const int numLocal = X.localNodes.size();
    localNodes.resize( numLocal );
    for( int s=0; s<numLocal; ++s )
        localNodes[s] = X.localNodes[s];

    return *this;
}

template<typename F>
inline void
DistNodalMatrix<F>::Pull
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::Pull");
#endif
    DistNodalMultiVec<F> XMultiVec( inverseMap, info, X );
    *this = XMultiVec;
    ComputeCommMetas( info );
}

template<typename F>
inline void
DistNodalMatrix<F>::Push
( const DistMap& inverseMap, const DistSymmInfo& info,
        DistMultiVec<F>& X ) const
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::Push");
#endif
    DistNodalMultiVec<F> XMultiVec( *this );
    XMultiVec.Push( inverseMap, info, X );
}

template<typename F>
inline int
DistNodalMatrix<F>::Height() const
{ return height_; }

template<typename F>
inline int
DistNodalMatrix<F>::Width() const
{ return width_; }

template<typename F>
inline void
DistNodalMatrix<F>::ComputeCommMetas( const DistSymmInfo& info ) const
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::ComputeCommMetas");
#endif
    const int numDist = info.distNodes.size();
    commMetas.resize( numDist-1 );

    // Handle the non-trivially distributed nodes
    for( int s=1; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        const int teamSize = mpi::CommSize( node.comm );
        const int teamRank = mpi::CommRank( node.comm );
        const Grid& grid = *node.grid;
        const int gridHeight = grid.Height();
        const int gridWidth = grid.Width();

        const DistSymmNodeInfo& childNode = info.distNodes[s-1];
        const int childTeamSize = mpi::CommSize( childNode.comm );
        const int childTeamRank = mpi::CommRank( childNode.comm );
        const bool inFirstTeam = ( childTeamRank == teamRank );
        const bool leftIsFirst = ( childNode.onLeft==inFirstTeam );
        const int leftTeamSize =
            ( childNode.onLeft ? childTeamSize : teamSize-childTeamSize );
        const int rightTeamSize = teamSize - leftTeamSize;
        const int leftTeamOffset = ( leftIsFirst ? 0 : rightTeamSize );
        const int rightTeamOffset = ( leftIsFirst ? leftTeamSize : 0 );
        const Grid& childGrid = *childNode.grid;
        const int childGridHeight = childGrid.Height();
        const int childGridWidth = childGrid.Width();

        // Fill numChildSendInd
        MatrixCommMeta& commMeta = commMetas[s-1];
        commMeta.Empty();
        commMeta.numChildSendInd.resize( teamSize );
        elem::MemZero( &commMeta.numChildSendInd[0], teamSize );
        const int updateSize = childNode.lowerStruct.size();
        const std::vector<int>& myRelInd =
            ( childNode.onLeft ? node.leftRelInd : node.rightRelInd );
        {
            const int colAlign = childNode.size % childGridHeight;
            const int colShift =
                Shift( childGrid.Row(), colAlign, childGridHeight );
            const int localHeight =
                Length( updateSize, colShift, childGridHeight );

            const int rowShift = childGrid.Col();
            const int localWidth = Length( width_, rowShift, childGridWidth );
            for( int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
            {
                const int iChild = colShift + iChildLoc*childGridHeight;
                const int destRow = myRelInd[iChild] % gridHeight;
                for( int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
                {
                    const int jChild = rowShift + jChildLoc*childGridWidth;
                    const int destCol = jChild % gridWidth;
                    const int destRank = destRow + destCol*gridHeight;
                    ++commMeta.numChildSendInd[destRank];
                }
            }
        }

        const int numLeftInd = node.leftRelInd.size();
        const int numRightInd = node.rightRelInd.size();
        std::vector<int> leftRowInd, rightRowInd;
        for( int i=0; i<numLeftInd; ++i )
            if( node.leftRelInd[i] % gridHeight == grid.Row() )
                leftRowInd.push_back( i );
        for( int i=0; i<numRightInd; ++i )
            if( node.rightRelInd[i] % gridHeight == grid.Row() )
                rightRowInd.push_back( i );

        // Get the child grid dimensions
        int childGridDims[4];
        GetChildGridDims( node, childNode, childGridDims );
        const int leftGridHeight = childGridDims[0];
        const int leftGridWidth = childGridDims[1];
        const int rightGridHeight = childGridDims[2];
        const int rightGridWidth = childGridDims[3];

        //
        // Compute the solve recv indices
        //
        commMeta.childRecvInd.resize( teamSize );
        for( int q=0; q<teamSize; ++q )
            commMeta.childRecvInd[q].clear();
        const int colShift = grid.Row();
        const int rowShift = grid.Col();
        const int localWidth = Length( width_, rowShift, gridWidth );
        // Append the indices from the left child
        const int numLeftRowInd = leftRowInd.size();
        for( int iPre=0; iPre<numLeftRowInd; ++iPre )
        {
            const int iChild = leftRowInd[iPre];
            const int iFront = node.leftRelInd[iChild];
            const int iFrontLoc = (iFront-colShift) / gridHeight;
            const int childRow = (node.leftSize+iChild) % leftGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*gridWidth;
                const int childCol = j % leftGridWidth;
                const int childRank = childRow + childCol*leftGridHeight;
                const int frontRank = leftTeamOffset + childRank;
                commMeta.childRecvInd[frontRank].push_back(iFrontLoc);
                commMeta.childRecvInd[frontRank].push_back(jLoc);
             }
        }
        // Append the indices from the right child
        const int numRightRowInd = rightRowInd.size();
        for( int iPre=0; iPre<numRightRowInd; ++iPre )
        {
            const int iChild = rightRowInd[iPre];
            const int iFront = node.rightRelInd[iChild];
            const int iFrontLoc = (iFront-colShift) / gridHeight;
            const int childRow = (node.rightSize+iChild) % rightGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*gridWidth;
                const int childCol = j % rightGridWidth;
                const int childRank = childRow + childCol*rightGridHeight;
                const int frontRank = rightTeamOffset + childRank;
                commMeta.childRecvInd[frontRank].push_back(iFrontLoc);
                commMeta.childRecvInd[frontRank].push_back(jLoc);
            }
        }
    }
}

} // namespace cliq
