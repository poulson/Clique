/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERSOLVE_LOCAL_HPP
#define CLIQ_NUMERIC_LOWERSOLVE_LOCAL_HPP

namespace cliq {

template<typename F> 
void LocalLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X );
template<typename F> 
void LocalLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X );

template<typename F> 
void LocalLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X,
  bool conjugate=false );
template<typename F> 
void LocalLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X,
  bool conjugate=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerForwardSolve"))
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D || 
                            L.frontType == BLOCK_LDL_INTRAPIV_2D );
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];
        elem::MakeZeros( WB );

        // Update using the children (if they exist)
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftInd = node.children[0];
            const int rightInd = node.children[1];
            Matrix<F>& leftWork = L.localFronts[leftInd].work;
            Matrix<F>& rightWork = L.localFronts[rightInd].work;
            const int leftNodeSize = info.localNodes[leftInd].size;
            const int rightNodeSize = info.localNodes[rightInd].size;
            const int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            auto leftUpdate = 
                LockedView( leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = node.leftRelInds[iChild]; 
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            auto rightUpdate =
                LockedView
                ( rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = node.rightRelInds[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Solve against this front
        if( blockLDL )
            FrontBlockLowerForwardSolve( frontL, W );
        else
            FrontLowerForwardSolve( frontL, W );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }
}

// This is an exact copy of the DistNodalMultiVec version...
template<typename F> 
inline void LocalLowerForwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerForwardSolve"))
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D ||
                            L.frontType == BLOCK_LDL_INTRAPIV_2D );
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];
        elem::MakeZeros( WB );

        // Update using the children (if they exist)
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftInd = node.children[0];
            const int rightInd = node.children[1];
            Matrix<F>& leftWork = L.localFronts[leftInd].work;
            Matrix<F>& rightWork = L.localFronts[rightInd].work;
            const int leftNodeSize = info.localNodes[leftInd].size;
            const int rightNodeSize = info.localNodes[rightInd].size;
            const int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            auto leftUpdate =
                LockedView( leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = node.leftRelInds[iChild]; 
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            auto rightUpdate =
                LockedView
                ( rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = node.rightRelInds[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0

        // Solve against this front
        if( blockLDL )
            FrontBlockLowerForwardSolve( frontL, W );
        else
            FrontLowerForwardSolve( frontL, W );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }
}

template<typename F> 
inline void LocalLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerBackwardSolve"))
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D ||
                            L.frontType == BLOCK_LDL_INTRAPIV_2D );
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];

        // Update using the parent
        const int parent = node.parent;
        DEBUG_ONLY(
            if( parent < 0 )
            {
                std::ostringstream msg;
                msg << "Parent index was negative: " << parent;
                LogicError( msg.str() );
            }
            if( parent >= numLocalNodes )  
            {
                std::ostringstream msg;
                msg << "Parent index was too large: " << parent << " >= "
                    << numLocalNodes;
                LogicError( msg.str() );
            }
        )
        Matrix<F>& parentWork = L.localFronts[parent].work;
        const SymmNodeInfo& parentNode = info.localNodes[parent];
        const int currentUpdateSize = WB.Height();
        const std::vector<int>& parentRelInds = 
          ( node.onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelInds[iCurrent];
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( node.onLeft )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work1d.Empty();
        }

        // Solve against this front
        if( blockLDL )
            FrontBlockLowerBackwardSolve( frontL, W, conjugate );
        else
            FrontLowerBackwardSolve( frontL, W, conjugate );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    L.distFronts[0].work1d.Empty();
    for( int s=0; s<numLocalNodes; ++s )
        L.localFronts[s].work.Empty();
}

template<typename F> 
inline void LocalLowerBackwardSolve
( const DistSymmInfo& info, 
  const DistSymmFrontTree<F>& L, DistNodalMatrix<F>& X,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerBackwardSolve"))
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D ||
                            L.frontType == BLOCK_LDL_INTRAPIV_2D );
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];

        // Update using the parent
        const int parent = node.parent;
        DEBUG_ONLY(
            if( parent < 0 )
            {
                std::ostringstream msg;
                msg << "Parent index was negative: " << parent;
                LogicError( msg.str() );
            }
            if( parent >= numLocalNodes )  
            {
                std::ostringstream msg;
                msg << "Parent index was too large: " << parent << " >= "
                    << numLocalNodes;
                LogicError( msg.str() );
            }
        )
        Matrix<F>& parentWork = L.localFronts[parent].work;
        const SymmNodeInfo& parentNode = info.localNodes[parent];
        const int currentUpdateSize = WB.Height();
        const std::vector<int>& parentRelInds = 
          ( node.onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelInds[iCurrent];
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( node.onLeft )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work2d.Empty();
        }

        // Solve against this front
        if( blockLDL )
            FrontBlockLowerBackwardSolve( frontL, W, conjugate );
        else
            FrontLowerBackwardSolve( frontL, W, conjugate );

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    L.distFronts[0].work2d.Empty();
    for( int s=0; s<numLocalNodes; ++s )
        L.localFronts[s].work.Empty();
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERSOLVE_LOCAL_HPP
