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
#ifndef CLIQUE_NESTED_DISSECTION_HPP
#define CLIQUE_NESTED_DISSECTION_HPP 1

#ifdef HAVE_PARMETIS

#include "parmetis.h"

namespace cliq {

// TODO
/*
void NestedDissection
( const DistGraph& graph, DistSeparatorTree& sepTree, DistSymmElimTree& eTree );*/

// TODO
/*
// For 1 process
int Bisect
( const Graph& graph, Graph& leftChild, Graph& rightChild, 
  std::vector<int>& map );
*/

int Bisect
( const DistGraph& graph, DistGraph& child, 
  std::vector<int>& localMap, bool& haveLeftChild );

void MapIndices
( const std::vector<int>& localMap, int blocksize,
  std::vector<int>& indices, mpi::Comm comm );

void InvertMap
( const std::vector<int>& localMap, int blocksize,
   std::vector<int>& localInverseMap, mpi::Comm comm );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline int Bisect
( const DistGraph& graph, DistGraph& child, 
  std::vector<int>& localMap, bool& haveLeftChild )
{
#ifndef RELEASE
    PushCallStack("Bisect");
#endif
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );
    if( commSize == 1 )
        throw std::logic_error
        ("This routine assumes at least two processes are used, "
         "otherwise one child will be lost");

    // Describe the source distribution
    const int blocksize = graph.Blocksize();
    std::vector<idx_t> vtxDist( commSize+1 );
    for( int i=0; i<commSize; ++i )
        vtxDist[i] = i*blocksize;
    vtxDist[commSize] = graph.NumSources();

    // ParMETIS assumes that there are no self-connections, so we must
    // manually remove them from our graph
    const int numLocalEdges = graph.NumLocalEdges();
    int numLocalSelfEdges = 0;
    for( int i=0; i<numLocalEdges; ++i )
        if( graph.Source(i) == graph.Target(i) )
            ++numLocalSelfEdges;

    // Fill our local connectivity (ignoring self edges)
    const int numLocalValidEdges = numLocalEdges - numLocalSelfEdges;
    const int numLocalSources = graph.NumLocalSources();
    const int firstLocalSource = graph.FirstLocalSource();
    std::vector<idx_t> xAdj( numLocalSources+1 );
    std::vector<idx_t> adjacency( numLocalValidEdges );
    int validCounter=0;
    int sourceOffset=0;
    int prevSource=firstLocalSource-1;
    for( int localEdge=0; localEdge<numLocalEdges; ++localEdge )
    {
        const int source = graph.Source( localEdge );
        const int target = graph.Target( localEdge );
#ifndef RELEASE
        if( source < prevSource )
            throw std::runtime_error("sources were not properly sorted");
#endif
        while( source != prevSource )
        {
            xAdj[sourceOffset++] = validCounter;
            ++prevSource;
        }
        if( source != target )
        {
            adjacency[validCounter] = target;
            ++validCounter;
        }
    }
#ifndef RELEASE
    if( sourceOffset != numLocalSources )
        throw std::logic_error("Mistake in xAdj computation");
#endif
    xAdj[numLocalSources] = numLocalValidEdges;

    // Create space for the result
    localMap.resize( numLocalSources );

    mpi::Barrier( comm );
    if( commRank == 0 )
        std::cout << "Starting CliqBisect..." << std::endl;

    idx_t numParSeps = 10;
    idx_t numSeqSeps = 5;
    real_t imbalance = 1.1;
    idx_t sizes[3];
    const int retval = CliqBisect
    ( &vtxDist[0], &xAdj[0], &adjacency[0], &numParSeps, &numSeqSeps, 
      &imbalance, NULL, &localMap[0], sizes, &comm );

    const int leftChildSize = sizes[0];
    const int rightChildSize = sizes[1];
    const int sepSize = sizes[2];

    if( commRank == 0 )
    {
        std::cout << "sizes: " 
                  << sizes[0] << ", " << sizes[1] << ", " << sizes[2] 
                  << std::endl;
    }

    // Build the child graph from the partitioned parent
    const int smallTeamSize = commSize/2;
    const int largeTeamSize = commSize - smallTeamSize;
    const bool inSmallTeam = ( commRank < smallTeamSize );
    const bool smallOnLeft = ( leftChildSize <= rightChildSize );
    const int leftTeamSize = ( smallOnLeft ? smallTeamSize : largeTeamSize );
    const int rightTeamSize = ( smallOnLeft ? largeTeamSize : smallTeamSize );
    const int leftTeamOffset = ( smallOnLeft ? 0 : smallTeamSize );
    const int rightTeamOffset = ( smallOnLeft ? smallTeamSize : 0 );
    const int leftTeamBlocksize = leftChildSize / leftTeamSize;
    const int rightTeamBlocksize = rightChildSize / rightTeamSize;
    const bool inLeftTeam = ( smallOnLeft == inSmallTeam );

    // Count how many rows we must send to each process 
    std::vector<int> rowSendSizes( commSize, 0 );
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = localMap[s];
        if( i < leftChildSize )
        {
            const int q = leftTeamOffset + i/leftTeamBlocksize;
            ++rowSendSizes[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = 
                rightTeamOffset + (i-leftChildSize)/rightTeamBlocksize;
            ++rowSendSizes[q];
        }
    }

    // Exchange the number of rows
    std::vector<int> rowRecvSizes( commSize );
    mpi::AllToAll
    ( &rowSendSizes[0], 1,
      &rowRecvSizes[0], 1, comm );

    // Prepare for the AllToAll to exchange the row indices and 
    // the number of column indices per row
    int numSendRows=0;
    std::vector<int> rowSendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        rowSendOffsets[q] = numSendRows;
        numSendRows += rowSendSizes[q];
    }
    int numRecvRows=0;
    std::vector<int> rowRecvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        rowRecvOffsets[q] = numRecvRows;
        numRecvRows += rowRecvSizes[q];
    }

    // Pack the row indices and how many column entries there will be per row
    std::vector<int> rowSendLengths( numSendRows );
    std::vector<int> rowSendIndices( numSendRows );
    std::vector<int> offsets = rowSendOffsets;
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = localMap[s];
        if( i < leftChildSize )
        {
            const int q = leftTeamOffset + i/leftTeamBlocksize;
            rowSendLengths[offsets[q]] = i;
            rowSendIndices[offsets[q]] = graph.NumConnections( s );
            ++offsets[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = 
                rightTeamOffset + (i-leftChildSize)/rightTeamBlocksize;
            rowSendLengths[offsets[q]] = i;
            rowSendIndices[offsets[q]] = graph.NumConnections( s );
            ++offsets[q];
        }
    }

    // Perform the row lengths exchange
    std::vector<int> rowRecvLengths( numRecvRows );
    mpi::AllToAll
    ( &rowSendLengths[0], &rowSendSizes[0], &rowSendOffsets[0],
      &rowRecvLengths[0], &rowRecvSizes[0], &rowRecvOffsets[0], comm );
    rowSendLengths.clear();

    // Perform the row indices exchange
    std::vector<int> rowRecvIndices( numRecvRows );
    mpi::AllToAll
    ( &rowSendIndices[0], &rowSendSizes[0], &rowSendOffsets[0],
      &rowRecvIndices[0], &rowRecvSizes[0], &rowRecvOffsets[0], comm );
    rowSendIndices.clear();
    rowSendSizes.clear();
    rowSendOffsets.clear();

    // Set up for sending the column indices
    int numSendIndices=0;
    std::vector<int> indexSendSizes( commSize, 0 );
    std::vector<int> indexSendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowSendSizes[q];
        const int offset = rowSendOffsets[q];
        for( int s=0; s<numRows; ++s )
            indexSendSizes[q] += rowSendIndices[offset+s];

        indexSendOffsets[q] = numSendIndices;
        numSendIndices += indexSendSizes[q];
    }
    int numRecvIndices=0;
    std::vector<int> indexRecvSizes( commSize, 0 );
    std::vector<int> indexRecvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowRecvSizes[q];
        const int offset = rowRecvOffsets[q];
        for( int s=0; s<numRows; ++s )
            indexRecvSizes[q] += rowRecvIndices[offset+s];

        indexRecvOffsets[q] = numRecvIndices;
        numRecvIndices += indexRecvSizes[q];
    }

    // Pack the indices
    std::vector<int> sendIndices( numSendIndices );
    offsets = indexSendOffsets;
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = localMap[s];
        if( i < leftChildSize )
        {
            const int q = leftTeamOffset + i/leftTeamBlocksize;

            int& offset = offsets[q];
            const int numConnections = graph.NumConnections( s );
            const int localEdgeOffset = graph.LocalEdgeOffset( s );
            for( int j=0; j<numConnections; ++j )
                sendIndices[offset++] = graph.Target( localEdgeOffset+j );
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q =
                rightTeamOffset + (i-leftChildSize)/rightTeamBlocksize;
               
            int& offset = offsets[q];
            const int numConnections = graph.NumConnections( s );
            const int localEdgeOffset = graph.LocalEdgeOffset( s );
            for( int j=0u; j<numConnections; ++j )
                sendIndices[offset++] = graph.Target( localEdgeOffset+j );
        }
    }

    // Send/recv the column indices
    std::vector<int> recvIndices( numRecvIndices );
    mpi::AllToAll
    ( &sendIndices[0], &indexSendSizes[0], &indexSendOffsets[0],
      &recvIndices[0], &indexRecvSizes[0], &indexRecvOffsets[0], comm );
    sendIndices.clear();
    indexSendSizes.clear();
    indexSendOffsets.clear();

    // Get the indices after reordering
    mpi::Barrier( comm );
    if( commRank == 0 )
        std::cout << "Starting MapIndices..." << std::endl;
    MapIndices( localMap, blocksize, recvIndices, comm );

    // Put the connections into our new graph
    mpi::Barrier( comm );
    if( commRank == 0 )
        std::cout << "Starting assembly..." << std::endl;

    const int childTeamRank = 
        ( inLeftTeam ? commRank-leftTeamOffset : commRank-rightTeamOffset );
    MPI_Comm childComm;
    mpi::CommSplit( comm, inLeftTeam, childTeamRank, childComm );
    child.SetComm( childComm );
    if( inLeftTeam )
    {
        haveLeftChild = true;
        child.ResizeTo( leftChildSize );
    }
    else
    {
        haveLeftChild = false;
        child.ResizeTo( rightChildSize );
    }
    const int childFirstLocalSource = child.FirstLocalSource();
    child.StartAssembly();
    child.Reserve( recvIndices.size() );
    int offset=0;
    for( int s=0; s<numRecvRows; ++s )
    {
        const int source = rowRecvIndices[s];
        const int numConnections = rowRecvLengths[s];
        for( int t=0; t<numConnections; ++t )
        {
            const int target = recvIndices[offset++];
            child.PushBack( source, target );
        }
    }
    child.StopAssembly();
#ifndef RELEASE
    PopCallStack();
#endif
    return sepSize;
}

// Overwrite the array of indices with the distributed map defined by each 
// processes's localMap. 
inline void MapIndices
( const std::vector<int>& localMap, int blocksize,
  std::vector<int>& indices, mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("MapIndices");
#endif
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    const int firstLocalSource = blocksize*commSize;
    const int numLocalSources = localMap.size();
    const int numLocalIndices = indices.size();

    // Count how many indices we need each process to map
    std::vector<int> requestSizes( commSize, 0 );
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = indices[s];
        const int q = i / blocksize;
        ++requestSizes[q];
    }

    // Send our requests and find out what we need to fulfill
    std::vector<int> fulfillSizes( commSize );
    mpi::AllToAll
    ( &requestSizes[0], 1, 
      &fulfillSizes[0], 1, comm );

    // Prepare for the AllToAll to exchange request sizes
    int numRequests=0;
    std::vector<int> requestOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        requestOffsets[q] = numRequests;
        numRequests += requestSizes[q];
    }
#ifndef RELEASE
    if( numRequests != numLocalIndices )
        throw std::logic_error("Miscalculated numRequests");
#endif
    int numFulfills=0;
    std::vector<int> fulfillOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        fulfillOffsets[q] = numFulfills;
        numFulfills += fulfillSizes[q];
    }

    // Pack the requested information 
    std::vector<int> requests( numRequests );
    std::vector<int> offsets = requestOffsets;
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = indices[s];
        const int q = i / blocksize;
        requests[offsets[q]++] = i;
    }

    // Perform the first index exchange
    std::vector<int> fulfills( numFulfills );
    mpi::AllToAll
    ( &requests[0], &requestSizes[0], &requestOffsets[0],
      &fulfills[0], &fulfillSizes[0], &fulfillOffsets[0], comm );

    // Map all of the indices in 'fulfills'
    for( int s=0; s<numFulfills; ++s )
    {
        const int i = fulfills[s];
        const int iLocal = i - firstLocalSource;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= numLocalSources )
            throw std::logic_error("Invalid request");
#endif
        fulfills[s] = localMap[iLocal];
    }

    // Send everything back
    mpi::AllToAll
    ( &fulfills[0], &fulfillSizes[0], &fulfillOffsets[0],
      &requests[0], &requestSizes[0], &requestOffsets[0], comm );

    // Unpack in the same way we originally packed
    offsets = requestOffsets;
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = indices[s];
        const int q = i / blocksize;
        indices[s] = requests[offsets[q]++];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Generate our local portion of the inverse of a distributed map
inline void InvertMap
( const std::vector<int>& localMap, int blocksize,
  std::vector<int>& localInverseMap, mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("InvertMap");
#endif
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    const int firstLocalSource = blocksize*commSize;
    const int numLocalSources = localMap.size();

    // How many pairs of original and mapped indices to send to each process
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = localMap[s];
        const int q = i / blocksize;
        sendSizes[q] += 2;
    }

    // Coordinate all of the processes on their send sizes
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll
    ( &sendSizes[0], 1,
      &recvSizes[0], 1, comm );

    // Prepare for the AllToAll to exchange send sizes
    int numSends=0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSends;
        numSends += sendSizes[q];
    }
#ifndef RELEASE
    if( numSends != 2*numLocalSources )
        throw std::logic_error("Miscalculated numSends");
#endif
    int numReceives=0;
    std::vector<int> recvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvOffsets[q] = numReceives;
        numReceives += recvSizes[q];
    }
#ifndef RELEASE
    if( numReceives != 2*numLocalSources )
        throw std::logic_error("Mistake in number of receives");
#endif

    // Pack our map information
    std::vector<int> sends( numSends );
    std::vector<int> offsets = sendOffsets;
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = localMap[s];
        const int q = i / blocksize;
        sends[offsets[q]++] = s+firstLocalSource;
        sends[offsets[q]++] = i;
    }

    // Send out the map information
    std::vector<int> recvs( numReceives );
    mpi::AllToAll
    ( &sends[0], &sendSizes[0], &sendOffsets[0],
      &recvs[0], &recvSizes[0], &recvOffsets[0], comm );

    // Form our part of the inverse map
    localInverseMap.resize( numLocalSources );
    for( int s=0; s<numReceives; s+=2 )
    {
        const int origIndex = recvs[s];
        const int mappedIndex = recvs[s+1];
        localInverseMap[mappedIndex-firstLocalSource] = origIndex;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // HAVE_PARMETIS

#endif // CLIQUE_NESTED_DISSECTION_HPP
