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
    const int numLocalSupernodes = S.local.supernodes.size();
    for( int s=0; s<numLocalSupernodes; ++s )
    {
        LocalSymmFactSupernode& sn = S.local.supernodes[s];
        Matrix<F>& frontL = L.local.fronts[s].frontL;
        Matrix<F>& frontR = L.local.fronts[s].frontR;
#ifndef RELEASE
        if( frontL.Height() != sn.size+sn.lowerStruct.size() ||
            frontL.Width() != sn.size )
            throw std::logic_error("Front was not the proper size");
#endif
        frontR.ResizeTo( frontL.Height(), sn.lowerStruct.size() );
        frontR.SetToZero();

        // Add updates from children (if they exist)
        const int numChildren = sn.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = sn.children[0];
            const int rightIndex = sn.children[1];
            Matrix<F>& leftFrontR = L.local.fronts[leftIndex].frontR;
            Matrix<F>& rightFrontR = L.local.fronts[rightIndex].frontR;
            const int leftSupernodeSize = S.local.supernodes[leftIndex].size;
            const int rightSupernodeSize = S.local.supernodes[rightIndex].size;
            const int leftUpdateSize = leftFrontR.Width();
            const int rightUpdateSize = rightFrontR.Width();
            
            Matrix<F> leftUpdate;
            leftUpdate.LockedView
            ( leftFrontR, 
              leftSupernodeSize, 0, leftUpdateSize, leftUpdateSize );

            Matrix<F> rightUpdate;
            rightUpdate.LockedView
            ( rightFrontR, 
              rightSupernodeSize, 0, rightUpdateSize, rightUpdateSize );

            // Add the left child's update matrix
            for( int jChild=0; jChild<leftUpdateSize; ++jChild )
            {
                const int jFront = sn.leftChildRelIndices[jChild];
                for( int iChild=0; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = sn.leftChildRelIndices[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
                    if( jFront < sn.size )
                        frontL.Update( iFront, jFront, value );
                    else
                        frontR.Update( iFront, jFront-sn.size, value );
                }
            }
            leftFrontR.Empty();

            // Add the right child's update matrix
            for( int jChild=0; jChild<rightUpdateSize; ++jChild )
            {
                const int jFront = sn.rightChildRelIndices[jChild];
                for( int iChild=0; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = sn.rightChildRelIndices[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
                    if( jFront < sn.size )
                        frontL.Update( iFront, jFront, value );
                    else
                        frontR.Update( iFront, jFront-sn.size, value );
                }
            }
            rightFrontR.Empty();
        }

        // Call the custom partial LDL
        LocalFrontLDL( orientation, frontL, frontR );
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

