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
#ifndef CLIQUE_LOCAL_LOWER_MULTIPLY
#define CLIQUE_LOCAL_LOWER_MULTIPLY 1

namespace cliq {

template<typename F>
void LocalLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const SymmInfo& info, const SymmFrontTree<F>& L, Matrix<F>& X );

template<typename F>
void LocalLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const SymmInfo& info, const SymmFrontTree<F>& L, Matrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset,
  const SymmInfo& info, const SymmFrontTree<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("LocalLowerMultiplyNormal");
#endif
    const int numLocalNodes = info.local.nodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const LocalSymmNodeInfo& node = info.local.nodes[s];
        const Matrix<F>& frontL = L.local.fronts[s].frontL;
        Matrix<F>& W = L.local.fronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
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
            Matrix<F>& leftWork = L.local.fronts[leftIndex].work;
            Matrix<F>& rightWork = L.local.fronts[rightIndex].work;
            const int leftNodeSize = info.local.nodes[leftIndex].size;
            const int rightNodeSize = info.local.nodes[rightIndex].size;
            const int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            Matrix<F> leftUpdate;
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
            Matrix<F> rightUpdate;
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

template<typename F> 
inline void LocalLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const SymmInfo& info, const SymmFrontTree<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("LocalLowerMultiplyTranspose");
#endif
    const int numLocalNodes = info.local.nodes.size();
    const int width = X.Width();

    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const LocalSymmNodeInfo& node = info.local.nodes[s];
        const Matrix<F>& frontL = L.local.fronts[s].frontL;
        Matrix<F>& W = L.local.fronts[s].work;

        // Set up a workspace
        W.ResizeTo( frontL.Height(), width );
        Matrix<F> WT, WB;
        elem::PartitionDown
        ( W, WT,
             WB, node.size );

        // Pull in the relevant information from the RHS
        Matrix<F> XT;
        XT.View( X, node.myOffset, 0, node.size, width );
        WT = XT;

        // Update using the parent's portion of the RHS
        const int parent = node.parent;
        const LocalSymmNodeInfo& parentNode = info.local.nodes[parent];
        Matrix<F>& parentWork = L.local.fronts[parent].work;
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
                L.dist.fronts[0].work1d.Empty();
        }

        // Make a copy of the unmodified RHS
        Matrix<F> XNode = W;

        // Multiply the (conjugate-)transpose of this block column of L against
        // this node's portion of the right-hand side.
        LocalFrontLowerMultiply( orientation, diag, diagOffset, frontL, XNode );

        // Store this node's portion of the result
        Matrix<F> XNodeT, XNodeB;
        elem::PartitionDown
        ( XNode, XNodeT,
                 XNodeB, node.size );
        XT = XNodeT;
        XNode.Empty();
    }
    L.dist.fronts[0].work1d.Empty();
    L.local.fronts.front().work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // CLIQUE_LOCAL_LOWER_MULTIPLY
