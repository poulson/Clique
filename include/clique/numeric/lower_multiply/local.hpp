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

template<typename T>
void LocalLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& X );

template<typename T>
void LocalLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T> 
inline void LocalLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("LocalLowerMultiplyNormal");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        const Matrix<T>& frontL = L.localFronts[s].frontL;
        Matrix<T>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<T> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<T> XT;
        XT.View( X, node.myOffset, 0, node.size, width );
        WT = XT;
        elem::MakeZeros( WB );

        // Multiply this block column of L against this node's portion of the
        // right-hand side and set W equal to the result
        LocalFrontLowerMultiply( NORMAL, diag, diagOffset, frontL, W );

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
            leftUpdate.LockedView
            ( leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = node.leftChildRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            Matrix<T> rightUpdate;
            rightUpdate.LockedView
            ( rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = node.rightChildRelIndices[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0 

        // Store this node's portion of the result
        XT = WT;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void LocalLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("LocalLowerMultiplyTranspose");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();

    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        const Matrix<T>& frontL = L.localFronts[s].frontL;
        Matrix<T>& W = L.localFronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<T> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<T> XT;
        XT.View( X, node.myOffset, 0, node.size, width );
        WT = XT;

        // Update using the parent's portion of the RHS
        const int parent = node.parent;
        const LocalSymmNodeInfo& parentNode = info.localNodes[parent];
        Matrix<T>& parentWork = L.localFronts[parent].work;
        const int currentUpdateSize = WB.Height();
        const std::vector<int>& parentRelIndices = 
            ( node.isLeftChild ? 
              parentNode.leftChildRelIndices :
              parentNode.rightChildRelIndices );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelIndices[iCurrent]; 
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so we can 
        // safely free the parent's work if this node is the left child
        if( node.isLeftChild )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work1d.Empty();
        }

        // Make a copy of the unmodified RHS
        Matrix<T> XNode = W;

        // Multiply the (conjugate-)transpose of this block column of L against
        // this node's portion of the right-hand side.
        LocalFrontLowerMultiply( orientation, diag, diagOffset, frontL, XNode );

        // Store this node's portion of the result
        Matrix<T> XNodeT, XNodeB;
        elem::PartitionDown
        ( XNode, XNodeT,
                 XNodeB, node.size );
        XT = XNodeT;
        XNode.Empty();
    }
    L.distFronts[0].work1d.Empty();
    L.localFronts.front().work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
