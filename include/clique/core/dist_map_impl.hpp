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

inline
DistMap::DistMap()
: numSources_(0), comm_(mpi::COMM_WORLD)
{ SetComm( mpi::COMM_WORLD ); } 

inline
DistMap::DistMap( mpi::Comm comm )
: numSources_(0), comm_(mpi::COMM_WORLD)
{ SetComm( comm ); }

inline
DistMap::DistMap( int numSources, mpi::Comm comm )
: numSources_(numSources), comm_(mpi::COMM_WORLD)
{ SetComm( comm ); }

inline
DistMap::~DistMap()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommFree( comm_ ); 
}

inline void
DistMap::StoreOwners
( int numSources, std::vector<int>& localIndices, mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("DistMap::StoreOwners");
#endif
    SetComm( comm );
    ResizeTo( numSources );
    const int commSize = mpi::CommSize( comm );
    const int blocksize = Blocksize();
    const int numLocalSources = NumLocalSources();
    const int firstLocalSource = FirstLocalSource();

    // Exchange via AllToAlls
    std::vector<int> sendSizes( commSize, 0 );
    const int numLocalIndices = localIndices.size();
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = localIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++sendSizes[q];
    }
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm );
    std::vector<int> sendOffsets( commSize ),
                     recvOffsets( commSize );
    int numSends=0, numRecvs=0;
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSends;
        recvOffsets[q] = numRecvs;
        numSends += sendSizes[q];
        numRecvs += recvSizes[q];
    }
#ifndef RELEASE
    if( numRecvs != numLocalSources )
        throw std::logic_error("Incorrect number of recv indices");
#endif
    std::vector<int> offsets = sendOffsets;
    std::vector<int> sendIndices( numSends );
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = localIndices[s];
        const int q = RowToProcess( i, blocksize, commSize );
        sendIndices[offsets[q]++] = i;
    }
    std::vector<int> recvIndices( numRecvs );
    mpi::AllToAll
    ( &sendIndices[0], &sendSizes[0], &sendOffsets[0],
      &recvIndices[0], &recvSizes[0], &recvOffsets[0], comm );

    // Form map
    for( int q=0; q<commSize; ++q )
    {
        const int size = recvSizes[q];
        const int offset = recvOffsets[q];
        for( int s=0; s<size; ++s )
        {
            const int i = recvIndices[offset+s];
            const int iLocal = i - firstLocalSource;
            SetLocal( iLocal, q );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistMap::Translate( std::vector<int>& localIndices ) const
{
#ifndef RELEASE
    PushCallStack("DistMap::Translate");
#endif
    const int commSize = mpi::CommSize( comm_ );
    const int numLocalIndices = localIndices.size();

    // Count how many indices we need each process to map
    std::vector<int> requestSizes( commSize, 0 );
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = localIndices[s];
#ifndef RELEASE
        if( i < 0 )
            throw std::logic_error("Index was negative");
#endif
        if( i < numSources_ )
        {
            const int q = RowToProcess( i, blocksize_, commSize );
            ++requestSizes[q];
        }
    }

    // Send our requests and find out what we need to fulfill
    std::vector<int> fulfillSizes( commSize );
    mpi::AllToAll( &requestSizes[0], 1, &fulfillSizes[0], 1, comm_ );

    // Prepare for the AllToAll to exchange request sizes
    int numRequests=0;
    std::vector<int> requestOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        requestOffsets[q] = numRequests;
        numRequests += requestSizes[q];
    }
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
        const int i = localIndices[s];
        if( i < numSources_ )
        {
            const int q = RowToProcess( i, blocksize_, commSize );
            requests[offsets[q]++] = i;
        }
    }

    // Perform the first index exchange
    std::vector<int> fulfills( numFulfills );
    mpi::AllToAll
    ( &requests[0], &requestSizes[0], &requestOffsets[0],
      &fulfills[0], &fulfillSizes[0], &fulfillOffsets[0], comm_ );

    // Map all of the indices in 'fulfills'
    for( int s=0; s<numFulfills; ++s )
    {
        const int i = fulfills[s];
        const int iLocal = i - firstLocalSource_;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= (int)localMap_.size() )
        {
            const int commRank = mpi::CommRank( comm_ );
            std::ostringstream msg;
            msg << "invalid request: i=" << i << ", iLocal=" << iLocal
                << ", commRank=" << commRank << ", blocksize=" << blocksize_;
            throw std::logic_error( msg.str().c_str() );
        }
#endif
        fulfills[s] = localMap_[iLocal];
    }

    // Send everything back
    mpi::AllToAll
    ( &fulfills[0], &fulfillSizes[0], &fulfillOffsets[0],
      &requests[0], &requestSizes[0], &requestOffsets[0], comm_ );

    // Unpack in the same way we originally packed
    offsets = requestOffsets;
    for( int s=0; s<numLocalIndices; ++s )
    {
        const int i = localIndices[s];
        if( i < numSources_ )
        {
            const int q = RowToProcess( i, blocksize_, commSize );
            localIndices[s] = requests[offsets[q]++];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistMap::FormInverse( DistMap& inverseMap ) const
{
#ifndef RELEASE
    PushCallStack("DistMap::FormInverse");
#endif
    const int commSize = mpi::CommSize( comm_ );
    const int numLocalSources = localMap_.size();

    // How many pairs of original and mapped indices to send to each process
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = localMap_[s];
        const int q = RowToProcess( i, blocksize_, commSize );
        sendSizes[q] += 2;
    }

    // Coordinate all of the processes on their send sizes
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm_ );

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
        const int i = localMap_[s];
        const int q = RowToProcess( i, blocksize_, commSize );
        sends[offsets[q]++] = s+firstLocalSource_;
        sends[offsets[q]++] = i;
    }

    // Send out the map information
    std::vector<int> recvs( numReceives );
    mpi::AllToAll
    ( &sends[0], &sendSizes[0], &sendOffsets[0],
      &recvs[0], &recvSizes[0], &recvOffsets[0], comm_ );

    // Form our part of the inverse map
    inverseMap.numSources_ = numSources_;
    inverseMap.SetComm( comm_ );
    for( int s=0; s<numReceives; s+=2 )
    {
        const int origIndex = recvs[s];
        const int mappedIndex = recvs[s+1];
        inverseMap.SetLocal( mappedIndex-firstLocalSource_, origIndex );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistMap::Extend( DistMap& firstMap ) const
{
#ifndef RELEASE
    PushCallStack("DistMap::Extend");
    // TODO: Ensure that the communicators are congruent and that the maps
    //       are compatible sizes.
#endif
    Translate( firstMap.localMap_ ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistMap::Extend( const DistMap& firstMap, DistMap& compositeMap ) const
{
#ifndef RELEASE
    PushCallStack("DistMap::Extend");

#endif
    compositeMap = firstMap;
    Extend( compositeMap );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline int
DistMap::NumSources() const
{ return numSources_; }

inline void
DistMap::SetComm( mpi::Comm comm )
{
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommFree( comm_ );

    if( comm != mpi::COMM_WORLD )
        mpi::CommDup( comm, comm_ );
    else
        comm_ = comm;

    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    blocksize_ = numSources_/commSize;
    firstLocalSource_ = blocksize_*commRank;
    const int numLocalSources =
        ( commRank<commSize-1 ?
          blocksize_ :
          numSources_ - (commSize-1)*blocksize_ );
    localMap_.resize( numLocalSources );
}

inline mpi::Comm
DistMap::Comm() const
{ return comm_; }

inline int
DistMap::Blocksize() const
{ return blocksize_; }

inline int
DistMap::FirstLocalSource() const
{ return firstLocalSource_; }

inline int
DistMap::NumLocalSources() const
{ return localMap_.size(); }

inline int
DistMap::GetLocal( int localSource ) const
{ 
#ifndef RELEASE
    PushCallStack("DistMap::GetLocal");
    if( localSource < 0 || localSource >= (int)localMap_.size() )
        throw std::logic_error("local source is out of bounds");
#endif
    const int target = localMap_[localSource];
#ifndef RELEASE
    PopCallStack();
#endif
    return target;
}

inline void
DistMap::SetLocal( int localSource, int target )
{ 
#ifndef RELEASE
    PushCallStack("DistMap::SetLocal");
    if( localSource < 0 || localSource >= (int)localMap_.size() )
        throw std::logic_error("local source is out of bounds");
#endif
    localMap_[localSource] = target; 
#ifndef RELEASE
    PopCallStack();
#endif
}

inline std::vector<int>&
DistMap::LocalMap()
{ return localMap_; }

inline const std::vector<int>&
DistMap::LocalMap() const
{ return localMap_; }

inline int* 
DistMap::LocalBuffer() 
{ return &localMap_[0]; }

inline const int*
DistMap::LocalBuffer() const
{ return &localMap_[0]; }

inline void
DistMap::Empty()
{
    numSources_ = 0;
    blocksize_ = 0;
    firstLocalSource_ = 0;
    localMap_.clear();
}

inline void
DistMap::ResizeTo( int numSources )
{
    const int commRank = mpi::CommRank( comm_ );
    const int commSize = mpi::CommSize( comm_ );
    numSources_ = numSources;
    blocksize_ = numSources/commSize;
    firstLocalSource_ = commRank*blocksize_;
    const int numLocalSources = 
        ( commRank<commSize-1 ? blocksize_
                              : numSources-blocksize_*(commSize-1) );
    localMap_.resize( numLocalSources );
}

inline const DistMap&
DistMap::operator=( const DistMap& map )
{
    numSources_ = map.numSources_;
    SetComm( map.comm_ );
    localMap_ = map.localMap_;
    return *this;
}

} // namespace cliq
