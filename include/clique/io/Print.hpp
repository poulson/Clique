/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_IO_PRINT_HPP
#define CLIQ_IO_PRINT_HPP

namespace cliq {

void Print
( const Graph& graph, std::string msg="Graph", std::ostream& os=std::cout );
void Print
( const DistGraph& graph, std::string msg="DistGraph", 
  std::ostream& os=std::cout );

template<typename T>
void Print
( const SparseMatrix<T>& A, std::string msg="SparseMatrix", 
  std::ostream& os=std::cout );
template<typename T>
void Print
( const DistSparseMatrix<T>& A, std::string msg="DistSparseMatrix", 
  std::ostream& os=std::cout );

void PrintLocal( const DistSymmInfo& info, std::ostream& os=std::cout );

//
// Implementation begins here
//

inline void
Print( const Graph& graph, std::string msg, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print [Graph]"))
    if( msg != "" )
        os << msg << std::endl;
    const int numEdges = graph.NumEdges();
    const int* srcBuf = graph.LockedSourceBuffer();
    const int* tgtBuf = graph.LockedTargetBuffer();
    for( int e=0; e<numEdges; ++e )
        os << srcBuf[e] << " " << tgtBuf[e] << "\n";
    os << std::endl;
}

inline void
Print( const DistGraph& graph, std::string msg, std::ostream& os ) 
{
    DEBUG_ONLY(CallStackEntry cse("Print [DistGraph]"))
    const mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    const int numLocalEdges = graph.NumLocalEdges();
    std::vector<int> edgeSizes(commSize), edgeOffsets(commSize);
    mpi::AllGather( &numLocalEdges, 1, &edgeSizes[0], 1, comm );
    int numEdges=0;
    for( int q=0; q<commSize; ++q )
    {
        edgeOffsets[q] = numEdges;
        numEdges += edgeSizes[q];
    }

    std::vector<int> sources, targets;
    if( commRank == 0 )
    {
        sources.resize( numEdges );
        targets.resize( numEdges );
    }
    mpi::Gather
    ( graph.LockedSourceBuffer(), numLocalEdges,
      &sources[0], &edgeSizes[0], &edgeOffsets[0], 0, comm );
    mpi::Gather
    ( graph.LockedTargetBuffer(), numLocalEdges,
      &targets[0], &edgeSizes[0], &edgeOffsets[0], 0, comm );

    if( commRank == 0 )
    {
        if( msg != "" )
            os << msg << std::endl;
        for( int e=0; e<numEdges; ++e )
            os << sources[e] << " " << targets[e] << "\n";
        os << std::endl;
    }
}

template<typename T>
inline void
Print( const SparseMatrix<T>& A, std::string msg, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print [SparseMatrix]"))
    if( msg != "" )
        os << msg << std::endl;
    const int numEntries = A.NumEntries();
    const int* srcBuf = A.LockedSourceBuffer();
    const int* tgtBuf = A.LockedTargetBuffer();
    const T* valBuf = A.LockedValueBuffer();
    for( int s=0; s<numEntries; ++s )
        os << srcBuf[s] << " " << tgtBuf[s] << " " << valBuf[s] << "\n";
    os << std::endl;
}

template<typename T>
inline void
Print( const DistSparseMatrix<T>& A, std::string msg, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print [DistSparseMatrix]"))
    const mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    const int numLocalEntries = A.NumLocalEntries();
    std::vector<int> nonzeroSizes(commSize), nonzeroOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &nonzeroSizes[0], 1, comm );
    int numNonzeros=0;
    for( int q=0; q<commSize; ++q )
    {
        nonzeroOffsets[q] = numNonzeros;
        numNonzeros += nonzeroSizes[q];
    }

    std::vector<int> sources, targets;
    std::vector<T> values;
    if( commRank == 0 )
    {
        sources.resize( numNonzeros );
        targets.resize( numNonzeros );
        values.resize( numNonzeros );
    }
    mpi::Gather
    ( A.LockedSourceBuffer(), numLocalEntries,
      &sources[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );
    mpi::Gather
    ( A.LockedTargetBuffer(), numLocalEntries,
      &targets[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );
    mpi::Gather
    ( A.LockedValueBuffer(), numLocalEntries,
      &values[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );

    if( commRank == 0 )
    {
        if( msg != "" )
            os << msg << std::endl;
        for( int s=0; s<numNonzeros; ++s )
            os << sources[s] << " " << targets[s] << " " << values[s] << "\n";
        os << std::endl;
    }
}

inline void 
PrintLocal( const DistSymmInfo& info, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("PrintLocal [DistSymmInfo]"))
    os << "Local nodes:" << std::endl;
    const int numLocal = info.localNodes.size();
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        os << " size=" << node.size << ", offset=" << node.off << "\n"; 
    }

    os << "Distributed nodes:" << std::endl;
    const int numDist = info.distNodes.size();
    for( int s=0; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        os << " size=" << node.size << ", offset=" << node.off << "\n";
    }
}

} // namespace cliq

#endif // ifndef CLIQ_IO_PRINT_HPP
