/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
#include "clique.hpp"
using namespace elemental;

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDL
( Orientation orientation, symbolic::SymmFact& S, numeric::SymmFrontTree<F>& L )
{
    using namespace clique::symbolic;
#ifndef RELEASE
    PushCallStack("numeric::LocalLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = S.local.supernodes.size();
    for( int s=0; s<numSupernodes; ++s )
    {
        LocalSymmFactSupernode& sn = S.local.supernodes[s];
        Matrix<F>& front = L.local.fronts[s].front;
#ifndef RELEASE
        if( front.Height() != sn.size+sn.lowerStruct.size() ||
            front.Width()  != sn.size+sn.lowerStruct.size() )
            throw std::logic_error("Front was not the proper size");
#endif

        // Add updates from children (if they exist)
        const int numChildren = sn.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = sn.children[0];
            const int rightIndex = sn.children[1];
            const Matrix<F>& leftFront = L.local.fronts[leftIndex].front;
            const Matrix<F>& rightFront = L.local.fronts[rightIndex].front;
            const int leftSupernodeSize = S.local.supernodes[leftIndex].size;
            const int rightSupernodeSize = S.local.supernodes[rightIndex].size;
            const int leftUpdateSize = leftFront.Height()-leftSupernodeSize;
            const int rightUpdateSize = rightFront.Height()-rightSupernodeSize;
            
            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftFront, leftSupernodeSize, leftSupernodeSize,
              leftUpdateSize, leftUpdateSize );

            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightFront, rightSupernodeSize, rightSupernodeSize,
              rightUpdateSize, rightUpdateSize );

            // Add the left child's update matrix
            for( int jChild=0; jChild<leftUpdateSize; ++jChild )
            {
                const int jFront = sn.leftChildRelIndices[jChild];
                for( int iChild=0; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = sn.leftChildRelIndices[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
                    front.Update( iFront, jFront, value );
                }
            }

            // Add the right child's update matrix
            for( int jChild=0; jChild<rightUpdateSize; ++jChild )
            {
                const int jFront = sn.rightChildRelIndices[jChild];
                for( int iChild=0; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = sn.rightChildRelIndices[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
                    front.Update( iFront, jFront, value );
                }
            }
        }

        // Call the custom partial LDL
        LocalFrontLDL( orientation, front, sn.size );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDL
( Orientation orientation, 
  symbolic::SymmFact& S, numeric::SymmFrontTree<float>& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<double>& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<std::complex<float> >& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::SymmFact& S, numeric::SymmFrontTree<std::complex<double> >& L );

