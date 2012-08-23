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

namespace cliq {

typedef unsigned char byte;
 
bool Initialized();
void Initialize( int& argc, char**& argv );
void Finalize();

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif

// Pull in some of Elemental's imported libraries
namespace blas = elem::blas;
namespace lapack = elem::lapack;
namespace mpi = elem::mpi;

// Pull in a number of useful enums from Elemental
using namespace elem::distribution_wrapper;
using namespace elem::left_or_right_wrapper;
using namespace elem::orientation_wrapper;
using namespace elem::unit_or_non_unit_wrapper;
using namespace elem::upper_or_lower_wrapper;

// For scalar operations
using elem::Base;
using elem::Complex;
using elem::Abs;
using elem::Sqrt;

// Pull in a few classes from Elemental
using elem::Matrix;
using elem::Grid;
using elem::DistMatrix;

// Pull in a few indexing routines
using elem::Shift;
using elem::LocalLength;

inline int
RowToProcess( int i, int blocksize, int commSize )
{
    if( blocksize > 0 )
        return std::min( i/blocksize, commSize-1 );
    else
        return commSize-1;
}

inline int
Find( const std::vector<int>& sortedIndices, int index )
{
#ifndef RELEASE
    PushCallStack("Find");
#endif
    std::vector<int>::const_iterator vecIt;
    vecIt = std::lower_bound
        ( sortedIndices.begin(), sortedIndices.end(), index );
#ifndef RELEASE
    if( vecIt == sortedIndices.end() )
        throw std::logic_error("Could not find index");
#endif
    const int indexOffset = vecIt - sortedIndices.begin();
#ifndef RELEASE
    PopCallStack();
#endif
    return indexOffset;
}
    
inline int
Find( const std::vector<int>& sortedIndices, int index, std::string msg )
{
#ifndef RELEASE
    PushCallStack("Find");
#endif
    std::vector<int>::const_iterator vecIt;
    vecIt = std::lower_bound
        ( sortedIndices.begin(), sortedIndices.end(), index );
#ifndef RELEASE
    if( vecIt == sortedIndices.end() )
        throw std::logic_error( msg.c_str() );
#endif
    const int indexOffset = vecIt - sortedIndices.begin();
#ifndef RELEASE
    PopCallStack();
#endif
    return indexOffset;
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
    {
        if( actualRecvCounts[proc] != recvCounts[proc] )
        {
            std::ostringstream msg;
            msg << "Expected recv count of " << recvCounts[proc]
                << " but recv'd " << actualRecvCounts[proc]
                << " from process " << proc << "\n";
            throw std::logic_error( msg.str().c_str() );
        }
    }
    actualRecvCounts.clear();
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
            ( &recvBuffer[displ], count, proc, 0, comm,
              requests[rCount++] );
    }
# ifdef BARRIER_IN_ALLTOALLV
    // This should help ensure that recvs are posted before the sends
    mpi::Barrier( comm );
# endif
    for( int proc=0; proc<commSize; ++proc )
    {
        int count = sendCounts[proc];
        int displ = sendDispls[proc];
        if( count != 0 )
            mpi::ISend
            ( &sendBuffer[displ], count, proc, 0, comm,
              requests[rCount++] );
    }
    mpi::WaitAll( numSends+numRecvs, &requests[0], &statuses[0] );
    statuses.clear();
    requests.clear();
#else
    mpi::AllToAll
    ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
      &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm );
#endif
}

} // namespace cliq
