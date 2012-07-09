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

void NestedDissection
( const DistGraph& graph, DistSymmElimTree& eTree );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline void NestedDissection
( const DistGraph& graph, DistSymmElimTree& eTree )
{
#ifndef RELEASE
    PushCallStack("NestedDissection");
#endif
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    // ParMETIS only generates a full nested-dissection tree for powers of 
    // two numbers of processes, so, for now, that is all we will support
    unsigned temp = commSize;
    unsigned log2CommSize=0;
    while( temp >>= 1 )
        ++log2CommSize;
    if( 1u<<log2CommSize != commSize )
        throw std::runtime_error
        ("NestedDissection currently requires a power of two # of processes");

    // Describe the source distribution
    const unsigned blocksize = graph.Blocksize();
    std::vector<idx_t> vtxDist( commSize+1 );
    for( unsigned i=0; i<commSize; ++i )
        vtxDist[i] = i*blocksize;
    vtxDist[commSize] = graph.NumSources();

    // ParMETIS assumes that there are no self-connections, so we must
    // manually remove them from our graph
    const unsigned numLocalEdges = graph.NumLocalEdges();
    unsigned numLocalSelfEdges = 0;
    for( unsigned i=0; i<numLocalEdges; ++i )
        if( graph.Source(i) == graph.Target(i) )
            ++numLocalSelfEdges;

    // Fill our local connectivity (ignoring self edges)
    const unsigned numLocalValidEdges = numLocalEdges - numLocalSelfEdges;
    const unsigned numLocalSources = graph.NumLocalSources();
    const unsigned firstLocalSource = graph.FirstLocalSource();
    std::vector<idx_t> xAdj( numLocalSources+1 );
    std::vector<idx_t> adjacency( numLocalValidEdges );
    unsigned validCounter=0;
    unsigned sourceOffset=0;
    unsigned prevSource=firstLocalSource;
    xAdj[sourceOffset++] = 0;
    for( unsigned i=0; i<numLocalEdges; ++i )
    {
        const unsigned source = graph.Source( i );
        const unsigned target = graph.Target( i );
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
    std::vector<idx_t> sizes( 2*commSize );

    if( commRank == 0 )
        std::cout << "Starting ParMETIS_V32_NodeND..." << std::endl;

    idx_t numFlag = 0;
    idx_t numParSeps = 5;
    idx_t numSeqSeps = 5;
    real_t imbalance = 1.1;
    const int retval = ParMETIS_V32_NodeND
    ( &vtxDist[0], &xAdj[0], &adjacency[0], NULL, &numFlag, NULL, NULL,
      &numParSeps, &numSeqSeps, &imbalance, NULL, NULL, &order[0], &sizes[0],
      &comm );

    if( commRank == 0 )
    {
        if( retval == METIS_OK )
            std::cout << "ParMETIS_V32_NodeND finished successfully" 
                      << std::endl;
        else
            std::cout << "ParMETIS_V32_NodeND failed" << std::endl;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif // HAVE_PARMETIS

#endif // CLIQUE_NESTED_DISSECTION_HPP
