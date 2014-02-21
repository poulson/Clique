/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_CORE_ENVIRONMENT_IMPL_HPP
#define CLIQ_CORE_ENVIRONMENT_IMPL_HPP

#include "clique/core/environment/decl.hpp"

namespace cliq {

// For getting the MPI argument instance (for internal usage)
inline void Args::HandleVersion( std::ostream& os ) const
{
    std::string version = "--version";
    char** arg = std::find( argv_, argv_+argc_, version );
    const bool foundVersion = ( arg != argv_+argc_ );
    if( foundVersion )
    {
        if( mpi::WorldRank() == 0 )
            PrintVersion();
        throw elem::ArgException();
    }
}

inline void Args::HandleBuild( std::ostream& os ) const
{
    std::string build = "--build";
    char** arg = std::find( argv_, argv_+argc_, build );
    const bool foundBuild = ( arg != argv_+argc_ );
    if( foundBuild )
    {
        if( mpi::WorldRank() == 0 )
        {
            PrintVersion();
            PrintConfig();
            PrintCCompilerInfo();
            PrintCxxCompilerInfo();
        }
        throw elem::ArgException();
    }
}

// For processing command-line arguments
template<typename T>
inline T
Input( std::string name, std::string desc )
{ return GetArgs().Input<T>( name, desc ); }

template<typename T>
inline T
Input( std::string name, std::string desc, T defaultVal )
{ return GetArgs().Input( name, desc, defaultVal ); }

inline void
ProcessInput()
{ GetArgs().Process(); }

inline void
PrintInputReport()
{ GetArgs().PrintReport(); }

template<typename T>
inline bool IsSorted( const std::vector<T>& x )
{
    const Int vecLength = x.size();
    for( Int i=1; i<vecLength; ++i )
    {
        if( x[i] < x[i-1] )
            return false;
    }
    return true;
}

// While is_strictly_sorted exists in Boost, it does not exist in the STL (yet)
template<typename T>
inline bool IsStrictlySorted( const std::vector<T>& x )
{
    const Int vecLength = x.size();
    for( Int i=1; i<vecLength; ++i )
    {
        if( x[i] <= x[i-1] )
            return false;
    }
    return true;
}

inline void Union
( std::vector<Int>& both, 
  const std::vector<Int>& first, const std::vector<Int>& second )
{
    both.resize( first.size()+second.size() );
    std::vector<Int>::iterator it = std::set_union
      ( first.begin(), first.end(), second.begin(), second.end(), 
        both.begin() );
    both.resize( Int(it-both.begin()) );
}

inline std::vector<Int>
Union( const std::vector<Int>& first, const std::vector<Int>& second )
{ 
    std::vector<Int> both;
    Union( both, first, second );
    return both;
}

inline void RelativeIndices
( std::vector<Int>& relInds, 
  const std::vector<Int>& sub, const std::vector<Int>& full )
{
    const Int numSub = sub.size();
    relInds.resize( numSub );
    std::vector<Int>::const_iterator it = full.begin();
    for( Int i=0; i<numSub; ++i )
    {
        const Int index = sub[i];
        it = std::lower_bound( it, full.end(), index );
        DEBUG_ONLY(
            if( it == full.end() )
                LogicError("Index was not found");
        )
        relInds[i] = Int(it-full.begin());
    }
}

inline std::vector<Int> RelativeIndices
( const std::vector<Int>& sub, const std::vector<Int>& full )
{
    std::vector<Int> relInds;
    RelativeIndices( relInds, sub, full );
    return relInds;
}

inline Int
RowToProcess( Int i, Int blocksize, Int commSize )
{
    if( blocksize > 0 )
        return std::min( i/blocksize, commSize-1 );
    else
        return commSize-1;
}

inline Int
Find
( const std::vector<Int>& sortedInds, Int index, std::string msg )
{
    DEBUG_ONLY(CallStackEntry cse("Find"))
    std::vector<Int>::const_iterator vecIt;
    vecIt = std::lower_bound( sortedInds.begin(), sortedInds.end(), index );
    DEBUG_ONLY(
        if( vecIt == sortedInds.end() )
            LogicError( msg );
    )
    return vecIt - sortedInds.begin();
}

inline void
VerifySendsAndRecvs
( const std::vector<int>& sendCounts,
  const std::vector<int>& recvCounts, mpi::Comm comm )
{
    const int commSize = mpi::CommSize( comm );
    std::vector<int> actualRecvCounts(commSize);
    mpi::AllToAll
    ( &sendCounts[0],       1,
      &actualRecvCounts[0], 1, comm );
    for( int proc=0; proc<commSize; ++proc )
        if( actualRecvCounts[proc] != recvCounts[proc] )
            LogicError
            ("Expected recv count of ",recvCounts[proc],
             " but recv'd ",actualRecvCounts[proc]," from process ",proc);
}

template<typename T>
inline void
SparseAllToAll
( const std::vector<T>& sendBuffer,
  const std::vector<int>& sendCounts, const std::vector<int>& sendDispls,
        std::vector<T>& recvBuffer,
  const std::vector<int>& recvCounts, const std::vector<int>& recvDispls,
        mpi::Comm comm )
{
#ifdef USE_CUSTOM_ALLTOALLV
    const int commSize = mpi::CommSize( comm );
    int numSends=0,numRecvs=0;
    for( int proc=0; proc<commSize; ++proc )
    {
        if( sendCounts[proc] != 0 )
            ++numSends;
        if( recvCounts[proc] != 0 )
            ++numRecvs;
    }
    std::vector<mpi::Status> statuses(numSends+numRecvs);
    std::vector<mpi::Request> requests(numSends+numRecvs);
    int rCount=0;
    for( int proc=0; proc<commSize; ++proc )
    {
        int count = recvCounts[proc];
        int displ = recvDispls[proc];
        if( count != 0 )
            mpi::IRecv
            ( &recvBuffer[displ], count, proc, comm, requests[rCount++] );
    }
#ifdef BARRIER_IN_ALLTOALLV
    // This should help ensure that recvs are posted before the sends
    mpi::Barrier( comm );
#endif
    for( int proc=0; proc<commSize; ++proc )
    {
        int count = sendCounts[proc];
        int displ = sendDispls[proc];
        if( count != 0 )
            mpi::ISend
            ( &sendBuffer[displ], count, proc, comm, requests[rCount++] );
    }
    mpi::WaitAll( numSends+numRecvs, &requests[0], &statuses[0] );
#else
    mpi::AllToAll
    ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
      &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm );
#endif
}

} // namespace cliq

#endif // ifndef CLIQ_CORE_ENVIRONMENT_IMPL_HPP
