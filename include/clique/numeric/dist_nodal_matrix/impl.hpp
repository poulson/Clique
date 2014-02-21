/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_DISTNODALMATRIX_IMPL_HPP
#define CLIQ_NUMERIC_DISTNODALMATRIX_IMPL_HPP

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
    DEBUG_ONLY(CallStackEntry cse("DistNodalMatrix::DistNodalMatrix"))
    Pull( inverseMap, info, X );
}

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix( const DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistNodalMatrix::DistNodalMatrix"))
    *this = X;
    commMetas.clear();
}

template<typename F>
inline const DistNodalMatrix<F>&
DistNodalMatrix<F>::operator=( const DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistNodalMatrix::operator="))
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
    DEBUG_ONLY(CallStackEntry cse("DistNodalMatrix::Pull"))
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
    DEBUG_ONLY(CallStackEntry cse("DistNodalMatrix::Push"))
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
    DEBUG_ONLY(CallStackEntry cse("DistNodalMatrix::ComputeCommMetas"))
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
        const int leftTeamOff = ( leftIsFirst ? 0 : rightTeamSize );
        const int rightTeamOff = ( leftIsFirst ? leftTeamSize : 0 );
        const Grid& childGrid = *childNode.grid;
        const int childGridHeight = childGrid.Height();
        const int childGridWidth = childGrid.Width();

        // Fill numChildSendInds
        MatrixCommMeta& commMeta = commMetas[s-1];
        commMeta.Empty();
        commMeta.numChildSendInds.resize( teamSize );
        elem::MemZero( &commMeta.numChildSendInds[0], teamSize );
        const int updateSize = childNode.lowerStruct.size();
        const std::vector<int>& myRelInds =
            ( childNode.onLeft ? node.leftRelInds : node.rightRelInds );
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
                const int destRow = myRelInds[iChild] % gridHeight;
                for( int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
                {
                    const int jChild = rowShift + jChildLoc*childGridWidth;
                    const int destCol = jChild % gridWidth;
                    const int destRank = destRow + destCol*gridHeight;
                    ++commMeta.numChildSendInds[destRank];
                }
            }
        }

        const int numLeftInds = node.leftRelInds.size();
        const int numRightInds = node.rightRelInds.size();
        std::vector<int> leftRowInds, rightRowInds;
        for( int i=0; i<numLeftInds; ++i )
            if( node.leftRelInds[i] % gridHeight == grid.Row() )
                leftRowInds.push_back( i );
        for( int i=0; i<numRightInds; ++i )
            if( node.rightRelInds[i] % gridHeight == grid.Row() )
                rightRowInds.push_back( i );

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
        commMeta.childRecvInds.resize( teamSize );
        for( int q=0; q<teamSize; ++q )
            commMeta.childRecvInds[q].clear();
        const int colShift = grid.Row();
        const int rowShift = grid.Col();
        const int localWidth = Length( width_, rowShift, gridWidth );
        // Append the indices from the left child
        const int numLeftRowInds = leftRowInds.size();
        for( int iPre=0; iPre<numLeftRowInds; ++iPre )
        {
            const int iChild = leftRowInds[iPre];
            const int iFront = node.leftRelInds[iChild];
            const int iFrontLoc = (iFront-colShift) / gridHeight;
            const int childRow = (node.leftSize+iChild) % leftGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*gridWidth;
                const int childCol = j % leftGridWidth;
                const int childRank = childRow + childCol*leftGridHeight;
                const int frontRank = leftTeamOff + childRank;
                commMeta.childRecvInds[frontRank].push_back(iFrontLoc);
                commMeta.childRecvInds[frontRank].push_back(jLoc);
             }
        }
        // Append the indices from the right child
        const int numRightRowInds = rightRowInds.size();
        for( int iPre=0; iPre<numRightRowInds; ++iPre )
        {
            const int iChild = rightRowInds[iPre];
            const int iFront = node.rightRelInds[iChild];
            const int iFrontLoc = (iFront-colShift) / gridHeight;
            const int childRow = (node.rightSize+iChild) % rightGridHeight;
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*gridWidth;
                const int childCol = j % rightGridWidth;
                const int childRank = childRow + childCol*rightGridHeight;
                const int frontRank = rightTeamOff + childRank;
                commMeta.childRecvInds[frontRank].push_back(iFrontLoc);
                commMeta.childRecvInds[frontRank].push_back(jLoc);
            }
        }
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_DISTNODALMATRIX_IMPL_HPP
