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
void clique::symbolic::DistSymmetricFactorization
( const LocalSymmStructure& symmStruct, LocalSymmFact& symmFact )
{
#ifndef RELEASE
    PushCallStack("symbolic::DistSymmetricFactorization");
#endif
    mpi::Comm comm = symmStruct.comm;
    const unsigned commRank = mpi::CommRank( comm );
    const unsigned commSize = mpi::CommSize( comm );
    const unsigned depth = symmStruct.sizes.size();

#ifndef RELEASE
    // Use the naive algorithm for finding floor(log2(commSize)) since it
    // will be an ignorable portion of the overhead and will be more easily
    // extended to more general integer types.
    unsigned temp = commSize;
    unsigned log2CommSize = 0;
    while( temp >>= 1 )
        ++log2CommSize;
    if( log2CommSize != depth )
        throw std::runtime_error("Invalid distributed tree depth");
    if( 1u<<log2CommSize != commSize )
        throw std::runtime_error
        ("Power-of-two number of procs currently required");
#endif

    symmFact.sizes.resize( depth );
    symmFact.offsets.resize( depth );
    symmFact.lowerStructs.resize( depth );
    symmFact.leftChildMaps.resize( depth );
    symmFact.rightChildMaps.resize( depth );

    if( depth == 0 )
        return;

    // The lowest level is still serial
    symmFact.sizes[0] = symmStruct.sizes[0];
    symmFact.offsets[0] = symmStruct.offsets[0];
    symmFact.lowerStructs[0] = symmStruct.lowerStructs[0];
    symmFact.leftChildMaps[0].resize(0); // not needed at this level
    symmFact.rightChildMaps[0].resize(0); // not needed at this level

    // Perform the parallel symbolic factorization
    std::vector<int> sendBuffer, recvBuffer;
    std::vector<int> childLowerStructUnion, lowerStructUnion, fullStruct,
                    superNodeIndices;
    std::vector<int>::iterator it;
    for( unsigned s=1; s<depth; ++s )
    {
        // Determine our partner based upon the bits of 'commRank'
        const unsigned partner = commRank ^ (1u << s);
        const bool onLeft = ( (commRank & (1u << s)) == 0 ); 
        const std::vector<int>& myChildLowerStruct = symmFact.lowerStructs[s-1];
        const std::vector<int>& origLowerStruct = symmStruct.lowerStructs[s];
        std::vector<int>& lowerStruct = symmFact.lowerStructs[s];
        const int superNodeSize = symmStruct.sizes[s];
        const int superNodeOffset = symmStruct.offsets[s];

        // SendRecv the message lengths
        const int myChildLowerStructSize = myChildLowerStruct.size();
        int theirChildLowerStructSize;
        mpi::SendRecv
        ( &myChildLowerStructSize, 1, partner, 0,
          &theirChildLowerStructSize, 1, partner, 0, comm );
        // Perform the exchange
        sendBuffer.resize( myChildLowerStructSize );
        recvBuffer.resize( theirChildLowerStructSize );
        std::memcpy
        ( &sendBuffer[0], &myChildLowerStruct[0], 
          myChildLowerStructSize*sizeof(int) );
        mpi::SendRecv
        ( &sendBuffer[0], myChildLowerStructSize, partner, 0,
          &recvBuffer[0], theirChildLowerStructSize, partner, 0, comm );
        
        // Union the two child lower structures
        childLowerStructUnion.resize
        ( myChildLowerStructSize+theirChildLowerStructSize );
        std::set_union
        ( sendBuffer.begin(), sendBuffer.end(),
          recvBuffer.begin(), recvBuffer.end(), childLowerStructUnion.begin() );

        // Union the lower structure of this supernode
        lowerStructUnion.resize
        ( childLowerStructUnion.size()+origLowerStruct.size() );
        std::set_union
        ( childLowerStructUnion.begin(), childLowerStructUnion.end(),
          origLowerStruct.begin(), origLowerStruct.end(), 
          lowerStructUnion.begin() );

        // Union again with the supernode indices
        std::vector<int> superNodeIndices( superNodeSize );
        for( int i=0; i<superNodeSize; ++i )
            superNodeIndices[i] = superNodeOffset + i;
        fullStruct.resize( superNodeSize + lowerStructUnion.size() );
        std::set_union
        ( superNodeIndices.begin(), superNodeIndices.end(),
          lowerStructUnion.begin(), lowerStructUnion.end(),
          fullStruct.begin() );

        // Construct the mappings from the local indices to the new global
        // indices.
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
        symmFact.leftChildMaps[s].resize( numLeftIndices );
        it = fullStruct.begin();
        for( int i=0; i<numLeftIndices; ++i )
        {
            it = std::lower_bound( it, fullStruct.end(), leftIndices[i] );
            symmFact.leftChildMaps[s][i] = *it;
        }
        symmFact.rightChildMaps[s].resize( numRightIndices );
        it = fullStruct.begin();
        for( int i=0; i<numRightIndices; ++i )
        {
            it = std::lower_bound( it, fullStruct.end(), rightIndices[i] );
            symmFact.rightChildMaps[s][i] = *it;
        }

        // Form lower structure of this node by removing the supernode indices
        lowerStruct.resize( fullStruct.size() );
        std::set_difference
        ( fullStruct.begin(), fullStruct.end(),
          superNodeIndices.begin(), superNodeIndices.end(),
          lowerStruct.begin() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

