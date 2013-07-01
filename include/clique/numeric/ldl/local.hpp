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
void 
LocalLDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void 
LocalLDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL )
{
#ifndef RELEASE
    CallStackEntry entry("LocalLDL");
#endif
    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
    {
        SymmNodeInfo& node = info.localNodes[s];
        const int updateSize = node.lowerStruct.size();
        Matrix<F>& frontL = L.localFronts[s].frontL;
        Matrix<F>& frontBR = L.localFronts[s].work;
        frontBR.Empty();
#ifndef RELEASE
        if( frontL.Height() != node.size+updateSize ||
            frontL.Width() != node.size )
            throw std::logic_error("Front was not the proper size");
#endif

        // Add updates from children (if they exist)
        Zeros( frontBR, updateSize, updateSize );
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = node.children[0];
            const int rightIndex = node.children[1];
            Matrix<F>& leftUpdate = L.localFronts[leftIndex].work;
            Matrix<F>& rightUpdate = L.localFronts[rightIndex].work;

            // Add the left child's update matrix
            const int leftUpdateSize = leftUpdate.Height();
            for( int jChild=0; jChild<leftUpdateSize; ++jChild )
            {
                const int jFront = node.leftRelIndices[jChild];
                for( int iChild=0; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = node.leftRelIndices[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
                    if( jFront < node.size )
                        frontL.Update( iFront, jFront, value );
                    else if( iFront >= node.size )
                        frontBR.Update
                        ( iFront-node.size, jFront-node.size, value );
                }
            }
            leftUpdate.Empty();

            // Add the right child's update matrix
            const int rightUpdateSize = rightUpdate.Height();
            for( int jChild=0; jChild<rightUpdateSize; ++jChild )
            {
                const int jFront = node.rightRelIndices[jChild];
                for( int iChild=0; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = node.rightRelIndices[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
                    if( jFront < node.size )
                        frontL.Update( iFront, jFront, value );
                    else if( iFront >= node.size )
                        frontBR.Update
                        ( iFront-node.size, jFront-node.size, value );
                }
            }
            rightUpdate.Empty();
        }

        // Call the custom partial LDL
        if( !blockLDL )
        {
            if( L.isHermitian )
                FrontLDL( ADJOINT, frontL, frontBR );
            else
                FrontLDL( TRANSPOSE, frontL, frontBR );
            frontL.GetDiagonal( L.localFronts[s].diag );
            elem::SetDiagonal( frontL, F(1) );
        }
        else
        {
            if( L.isHermitian )
                FrontBlockLDL( ADJOINT, frontL, frontBR );
            else
                FrontBlockLDL( TRANSPOSE, frontL, frontBR );
        }
    }
}

} // namespace cliq
