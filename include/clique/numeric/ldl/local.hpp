/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename F> 
void LocalLDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F> 
inline void LocalLDL
( DistSymmInfo& info, DistSymmFrontTree<F>& L, bool blockLDL )
{
#ifndef RELEASE
    PushCallStack("LocalLDL");
#endif
    const int numLocalNodes = info.localNodes.size();
    for( int s=0; s<numLocalNodes; ++s )
    {
        LocalSymmNodeInfo& node = info.localNodes[s];
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
        elem::Zeros( updateSize, updateSize, frontBR );
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
                const int jFront = node.leftChildRelIndices[jChild];
                for( int iChild=0; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = node.leftChildRelIndices[iChild];
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
                const int jFront = node.rightChildRelIndices[jChild];
                for( int iChild=0; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = node.rightChildRelIndices[iChild];
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
                LocalFrontLDL( ADJOINT, frontL, frontBR );
            else
                LocalFrontLDL( TRANSPOSE, frontL, frontBR );
        }
        else
        {
            if( L.isHermitian )
                LocalFrontBlockLDL( ADJOINT, frontL, frontBR );
            else
                LocalFrontBlockLDL( TRANSPOSE, frontL, frontBR );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
