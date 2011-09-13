/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
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

// This is the part of the symbolic factorization that requires fine-grain 
// parallelism: we assume that the upper floor(log2(commSize)) levels of the
// tree are balanced.
//
// TODO: Generalize so that the depth can be less than or equal to 
// floor(log2(commSize)). This would allow for the use of more processes in the
// factorization.
//
// TODO: Generalize to support more than just a power-of-two number of 
//       processes. This should be relatively straightforward.
void clique::symbolic::TopSymmetricFactorization
( const TopOrigStruct& topOrig,
  const BottomFactStruct& bottomFact,
        TopFactStruct& topFact )
{
#ifndef RELEASE
    PushCallStack("symbolic::TopSymmetricFactorization");
#endif
    const unsigned numSupernodes = topOrig.lowerStructs.size();
    topFact.comm = topOrig.comm;
    topFact.sizes = topOrig.sizes;
    topFact.offsets = topOrig.offsets;
    topFact.lowerStructs.resize( numSupernodes );
    topFact.origLowerRelIndices.resize( numSupernodes );
    topFact.leftChildRelIndices.resize( numSupernodes );
    topFact.rightChildRelIndices.resize( numSupernodes );
    if( numSupernodes == 0 )
        return;

    const unsigned commRank = mpi::CommRank( topOrig.comm );
    const unsigned commSize = mpi::CommSize( topOrig.comm );
#ifndef RELEASE
    // Use the naive algorithm for finding floor(log2(commSize)) since it
    // will be an ignorable portion of the overhead and will be more easily
    // extended to more general integer types.
    unsigned temp = commSize;
    unsigned log2CommSize = 0;
    while( temp >>= 1 )
        ++log2CommSize;
    if( log2CommSize != numSupernodes )
        throw std::runtime_error("Invalid distributed tree depth");
    if( 1u<<log2CommSize != commSize )
        throw std::runtime_error
        ("Power-of-two number of procs currently required");
#endif

    // The bottom node is already computed, so just copy it over
    topFact.lowerStructs[0] = bottomFact.lowerStructs.back();
    topFact.origLowerRelIndices[0] = bottomFact.origLowerRelIndices.back();
    topFact.leftChildRelIndices[0] = bottomFact.leftChildRelIndices.back();
    topFact.rightChildRelIndices[0] = bottomFact.rightChildRelIndices.back();

    // Perform the parallel symbolic factorization
    std::vector<int> sendBuffer, recvBuffer;
    std::vector<int> childrenStruct, partialStruct, fullStruct,
                    supernodeIndices;
    for( unsigned k=1; k<numSupernodes; ++k )
    {
        const int supernodeSize = topOrig.sizes[k];
        const int supernodeOffset = topOrig.offsets[k];
        const std::vector<int>& origLowerStruct = topOrig.lowerStructs[k];
        const std::vector<int>& myChildLowerStruct = topFact.lowerStructs[k-1];
        std::vector<int>& lowerStruct = topFact.lowerStructs[k];

        // Determine our partner based upon the bits of 'commRank'
        const unsigned partner = commRank ^ (1u << (k-1));
        const bool onLeft = ( (commRank & (1u << (k-1))) == 0 ); 

        // SendRecv the message lengths
        const int myChildLowerStructSize = myChildLowerStruct.size();
        int theirChildLowerStructSize;
        mpi::SendRecv
        ( &myChildLowerStructSize, 1, partner, 0,
          &theirChildLowerStructSize, 1, partner, 0, topOrig.comm );
        // Perform the exchange
        sendBuffer.resize( myChildLowerStructSize );
        recvBuffer.resize( theirChildLowerStructSize );
        std::memcpy
        ( &sendBuffer[0], &myChildLowerStruct[0], 
          myChildLowerStructSize*sizeof(int) );
        mpi::SendRecv
        ( &sendBuffer[0], myChildLowerStructSize, partner, 0,
          &recvBuffer[0], theirChildLowerStructSize, partner, 0, topOrig.comm );
        
        // Union the two child lower structures
        childrenStruct.resize
        ( myChildLowerStructSize+theirChildLowerStructSize );
        std::set_union
        ( sendBuffer.begin(), sendBuffer.end(),
          recvBuffer.begin(), recvBuffer.end(), childrenStruct.begin() );

        // Union the lower structure of this supernode
        partialStruct.resize( childrenStruct.size() + origLowerStruct.size() );
        std::set_union
        ( childrenStruct.begin(), childrenStruct.end(),
          origLowerStruct.begin(), origLowerStruct.end(), 
          partialStruct.begin() );

        // Union again with the supernode indices
        supernodeIndices.resize( supernodeSize );
        for( int i=0; i<supernodeSize; ++i )
            supernodeIndices[i] = supernodeOffset + i;
        fullStruct.resize( supernodeSize + partialStruct.size() );
        std::set_union
        ( supernodeIndices.begin(), supernodeIndices.end(),
          partialStruct.begin(), partialStruct.end(), 
          fullStruct.begin() );

        // Construct the relative indices of the original lower structure
        const int numOrigLowerIndices = origLowerStruct.size();
        topFact.origLowerRelIndices[k].resize( numOrigLowerIndices );
        std::vector<int>::iterator it = fullStruct.begin();
        for( int i=0; i<numOrigLowerIndices; ++i )
        {
            it = std::lower_bound( it, fullStruct.end(), origLowerStruct[i] );
            topFact.origLowerRelIndices[k][i] = *it;
        }

        // Construct the relative indices of the children
        int numLeftIndices, numRightIndices;
        const int *leftIndices, *rightIndices;
        if( onLeft )
        {
            leftIndices = &sendBuffer[0];
            rightIndices = &recvBuffer[0];
            numLeftIndices = sendBuffer.size();
            numRightIndices = recvBuffer.size();
        }
        else
        {
            leftIndices = &recvBuffer[0];
            rightIndices = &sendBuffer[0];
            numLeftIndices = recvBuffer.size();
            numRightIndices = sendBuffer.size();
        }
        topFact.leftChildRelIndices[k].resize( numLeftIndices );
        topFact.rightChildRelIndices[k].resize( numRightIndices );
        it = fullStruct.begin();
        for( int i=0; i<numLeftIndices; ++i )
        {
            it = std::lower_bound( it, fullStruct.end(), leftIndices[i] );
            topFact.leftChildRelIndices[k][i] = *it;
        }
        it = fullStruct.begin();
        for( int i=0; i<numRightIndices; ++i )
        {
            it = std::lower_bound( it, fullStruct.end(), rightIndices[i] );
            topFact.rightChildRelIndices[k][i] = *it;
        }

        // Form lower structure of this node by removing the supernode indices
        lowerStruct.resize( fullStruct.size() );
        std::set_difference
        ( fullStruct.begin(), fullStruct.end(),
          supernodeIndices.begin(), supernodeIndices.end(),
          lowerStruct.begin() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

