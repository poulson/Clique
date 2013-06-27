/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename T>
void LocalLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X );

template<typename T>
void LocalLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, 
  DistNodalMultiVec<T>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T> 
inline void LocalLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("LocalLowerMultiplyNormal");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<T>& frontL = L.localFronts[s].frontL;
        Matrix<T>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<T> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, node.size );
        WT = X.localNodes[s];
        elem::MakeZeros( WB );

        // Multiply this block column of L against this node's portion of the
        // right-hand side and set W equal to the result
        FrontLowerMultiply( NORMAL, diag, diagOffset, frontL, W );

        // Update using the children (if they exist)
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = node.children[0];
            const int rightIndex = node.children[1];
            Matrix<T>& leftWork = L.localFronts[leftIndex].work;
            Matrix<T>& rightWork = L.localFronts[rightIndex].work;
            const int leftNodeSize = info.localNodes[leftIndex].size;
            const int rightNodeSize = info.localNodes[rightIndex].size;
            const int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            Matrix<T> leftUpdate;
            LockedView
            ( leftUpdate, leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = node.leftRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            Matrix<T> rightUpdate;
            LockedView
            ( rightUpdate, 
              rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = node.rightRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0 

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }
}

template<typename T> 
inline void LocalLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, 
  DistNodalMultiVec<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("LocalLowerMultiplyTranspose");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();

    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<T>& frontL = L.localFronts[s].frontL;
        Matrix<T>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<T> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, node.size );
        WT = X.localNodes[s];

        // Update using the parent's portion of the RHS
        const int parent = node.parent;
        const SymmNodeInfo& parentNode = info.localNodes[parent];
        Matrix<T>& parentWork = L.localFronts[parent].work;
        const int currentUpdateSize = WB.Height();
        const std::vector<int>& parentRelIndices = 
            ( node.onLeft ? parentNode.leftRelIndices
                          : parentNode.rightRelIndices );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelIndices[iCurrent]; 
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so we can 
        // safely free the parent's work if this node is the left child
        if( node.onLeft )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work1d.Empty();
        }

        // Make a copy of the unmodified RHS
        Matrix<T> XNode = W;

        // Multiply the (conjugate-)transpose of this block column of L against
        // this node's portion of the right-hand side.
        FrontLowerMultiply( orientation, diag, diagOffset, frontL, XNode );

        // Store this node's portion of the result
        Matrix<T> XNodeT, XNodeB;
        elem::PartitionDown
        ( XNode, XNodeT,
                 XNodeB, node.size );
        X.localNodes[s] = XNodeT;
        XNode.Empty();
    }
    L.distFronts[0].work1d.Empty();
    L.localFronts.front().work.Empty();
}

} // namespace cliq
