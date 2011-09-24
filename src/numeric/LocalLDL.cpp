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
( Orientation orientation,
  symbolic::LocalSymmFact& S, // can't be const due to map...
  numeric::LocalSymmFact<F>& L )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = S.supernodes.size();
    for( int s=0; s<numSupernodes; ++s )
    {
        symbolic::LocalSymmFactSupernode& symbSN = S.supernodes[s];
        numeric::LocalSymmFactSupernode<F>& numSN = L.supernodes[s];
#ifndef RELEASE
        if( numSN.front.Height() != symbSN.size+symbSN.lowerStruct.size() ||
            numSN.front.Width()  != symbSN.size+symbSN.lowerStruct.size() )
            throw std::logic_error("Front was not the proper size");
#endif

        // Add updates from children (if they exist)
        const int numChildren = symbSN.children.size();
        if( numChildren == 2 )
        {
            const int leftIndex = symbSN.children[0];
            const int rightIndex = symbSN.children[1];
            const Matrix<F>& leftFront = L.supernodes[leftIndex].front;
            const Matrix<F>& rightFront = L.supernodes[rightIndex].front;
            const int leftSupernodeSize = S.supernodes[leftIndex].size;
            const int rightSupernodeSize = S.supernodes[rightIndex].size;
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
                const int jFront = symbSN.leftChildRelIndices[jChild];
                for( int iChild=0; iChild<leftUpdateSize; ++iChild )
                {
                    const int iFront = symbSN.leftChildRelIndices[iChild];
                    const F value = leftUpdate.Get(iChild,jChild);
                    numSN.front.Update( iFront, jFront, value );
                }
            }

            // Add the right child's update matrix
            for( int jChild=0; jChild<rightUpdateSize; ++jChild )
            {
                const int jFront = symbSN.rightChildRelIndices[jChild];
                for( int iChild=0; iChild<rightUpdateSize; ++iChild )
                {
                    const int iFront = symbSN.rightChildRelIndices[iChild];
                    const F value = rightUpdate.Get(iChild,jChild);
                    numSN.front.Update( iFront, jFront, value );
                }
            }
        }

        // Call the custom partial LDL
        LocalSupernodeLDL( orientation, numSN.front, symbSN.size );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalSymmFact& S,
  numeric::LocalSymmFact<float>& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalSymmFact& S,
  numeric::LocalSymmFact<double>& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalSymmFact& S,
  numeric::LocalSymmFact<std::complex<float> >& L );

template void clique::numeric::LocalLDL
( Orientation orientation,
  symbolic::LocalSymmFact& S,
  numeric::LocalSymmFact<std::complex<double> >& L );

