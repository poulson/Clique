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

template<typename F> 
void LocalLowerForwardSolve
( UnitOrNonUnit diag, 
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& X );

template<typename F> 
void LocalLowerBackwardSolve
( Orientation orientation, UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalLowerForwardSolve
( UnitOrNonUnit diag, 
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("LocalLowerForwardSolve");
#endif
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D );
#ifndef RELEASE
    if( blockLDL && diag == UNIT )
        throw std::logic_error("Unit diagonal is nonsensical for block LDL");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& W = L.localFronts[s].work;

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

        // Update using the children (if they exist)
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = node.children[0];
            const int rightIndex = node.children[1];
            Matrix<F>& leftWork = L.localFronts[leftIndex].work;
            Matrix<F>& rightWork = L.localFronts[rightIndex].work;
            const int leftNodeSize = info.localNodes[leftIndex].size;
            const int rightNodeSize = info.localNodes[rightIndex].size;
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

        // Solve against this front
        if( !blockLDL )
            LocalFrontLowerForwardSolve( diag, frontL, W );
        else
            LocalFrontBlockLowerForwardSolve( frontL, W );

        // Store this node's portion of the result
        XT = WT;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void LocalLowerBackwardSolve
( Orientation orientation, UnitOrNonUnit diag,
  const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("LocalLowerBackwardSolve");
#endif
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D );
#ifndef RELEASE
    if( blockLDL && diag == UNIT )
        throw std::logic_error("Unit diagonal is nonsensical for block LDL");
#endif
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const LocalSymmNodeInfo& node = info.localNodes[s];
        const Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& W = L.localFronts[s].work;

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

        // Update using the parent
        const int parent = node.parent;
#ifndef RELEASE
        if( parent < 0 )
        {
            std::ostringstream msg;
            msg << "Parent index was negative: " << parent;
            throw std::logic_error( msg.str().c_str() );
        }
        if( parent >= numLocalNodes )  
        {
            std::ostringstream msg;
            msg << "Parent index was too large: " << parent << " >= "
                << numLocalNodes;
            throw std::logic_error( msg.str().c_str() );
        }
#endif
        Matrix<F>& parentWork = L.localFronts[parent].work;
        const LocalSymmNodeInfo& parentNode = info.localNodes[parent];
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

        // The left child is numbered lower than the right child, so 
        // we can safely free the parent's work if we are the left child
        if( node.isLeftChild )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work1d.Empty();
        }

        // Solve against this front
        if( !blockLDL )
            LocalFrontLowerBackwardSolve( orientation, diag, frontL, W );
        else
            LocalFrontBlockLowerBackwardSolve( orientation, frontL, W );

        // Store this node's portion of the result
        XT = WT;
    }

    // Ensure that all of the temporary buffers are freed (this is overkill)
    L.distFronts[0].work1d.Empty();
    for( int s=0; s<numLocalNodes; ++s )
        L.localFronts[s].work.Empty();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
