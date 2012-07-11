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

// NOTE: This routine is nowhere near finished
void NestedDissection
( const DistGraph& graph, DistSeparatorTree& sepTree, DistSymmElimTree& eTree );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline void NestedDissection
( const DistGraph& graph, DistSeparatorTree& sepTree, DistSymmElimTree& eTree )
{
#ifndef RELEASE
    PushCallStack("NestedDissection");
#endif
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

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
    std::vector<idx_t> order( numLocalSources );

    if( commRank == 0 )
        std::cout << "Starting CliqBisect..." << std::endl;

    idx_t numParSeps = 5;
    idx_t numSeqSeps = 5;
    real_t imbalance = 1.1;
    idx_t size;
    const int retval = CliqBisect
    ( &vtxDist[0], &xAdj[0], &adjacency[0], &numParSeps, &numSeqSeps, 
      &imbalance, NULL, &order[0], &size, &comm );

    if( commRank == 0 )
    {
        if( retval == METIS_OK )
            std::cout << "CliqBisect was successful with size=" << size
                      << std::endl;
        else
            std::cout << "CliqBisect failed" << std::endl;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // HAVE_PARMETIS

#endif // CLIQUE_NESTED_DISSECTION_HPP
