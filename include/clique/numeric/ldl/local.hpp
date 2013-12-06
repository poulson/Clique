/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LDL_LOCAL_HPP
#define CLIQ_NUMERIC_LDL_LOCAL_HPP

namespace cliq {

template<typename F> 
void 
LocalLDL( DistSymmInfo& info, DistSymmFrontTree<F>& L );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void 
LocalLDL( DistSymmInfo& info, DistSymmFrontTree<F>& L )
{
#ifndef RELEASE
    CallStackEntry cse("LocalLDL");
#endif
    const bool blockLDL = ( L.frontType == BLOCK_LDL_2D ||
                            L.frontType == BLOCK_LDL_INTRAPIV_2D );
    const bool intraPiv = ( L.frontType == BLOCK_LDL_INTRAPIV_2D );

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
            LogicError("Front was not the proper size");
#endif

        // Add updates from children (if they exist)
        Zeros( frontBR, updateSize, updateSize );
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftInd = node.children[0];
            const int rightInd = node.children[1];
            Matrix<F>& leftUpdate = L.localFronts[leftInd].work;
            Matrix<F>& rightUpdate = L.localFronts[rightInd].work;

            // Add the left child's update matrix
            const int leftUpdateSize = leftUpdate.Height();
            for( int jChild=0; jChild<leftUpdateSize; ++jChild )
            {
                const int jFront = node.leftRelInds[jChild];
                for( int iChild=jChild; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = node.leftRelInds[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
#ifndef RELEASE
                    if( iFront < jFront )
                        LogicError("Tried to update upper triangle");
#endif
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
                const int jFront = node.rightRelInds[jChild];
                for( int iChild=jChild; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = node.rightRelInds[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
#ifndef RELEASE
                    if( iFront < jFront )
                        LogicError("Tried to update upper triangle");
#endif
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
        if( blockLDL )
            FrontBlockLDL( frontL, frontBR, L.isHermitian, intraPiv );
        else
        {
            FrontLDL( frontL, frontBR, L.isHermitian );
            frontL.GetDiagonal( L.localFronts[s].diag );
            elem::SetDiagonal( frontL, F(1) );
        }
    }
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LDL_LOCAL_HPP
