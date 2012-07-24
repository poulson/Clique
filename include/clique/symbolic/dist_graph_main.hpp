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
#ifndef CLIQUE_DIST_GRAPH_MAIN_HPP
#define CLIQUE_DIST_GRAPH_MAIN_HPP 1

namespace cliq {

inline 
DistGraph::DistGraph()
: numSources_(0), numTargets_(0)
{
    SetComm( mpi::COMM_WORLD );
}

inline
DistGraph::DistGraph( mpi::Comm comm )
: numSources_(0), numTargets_(0)
{
    SetComm( comm );
}

inline 
DistGraph::DistGraph( int numVertices, mpi::Comm comm )
: numSources_(numVertices), numTargets_(numVertices)
{
    SetComm( comm );
}

inline 
DistGraph::DistGraph( int numSources, int numTargets, mpi::Comm comm )
: numSources_(numSources), numTargets_(numTargets)
{
    SetComm( comm );
}

inline
DistGraph::DistGraph( const Graph& graph )
{
#ifndef RELEASE
    PushCallStack("DistGraph::DistGraph");
#endif
    *this = graph;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline
DistGraph::DistGraph( const DistGraph& graph )
{
#ifndef RELEASE
    PushCallStack("DistGraph::DistGraph");
#endif
    if( &graph != this )
        *this = graph;
    else
        throw std::logic_error("Tried to construct DistGraph with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

inline
DistGraph::~DistGraph()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommFree( comm_ );
} 

inline int 
DistGraph::NumSources() const
{ return numSources_; }

inline int 
DistGraph::NumTargets() const
{ return numTargets_; }

inline void 
DistGraph::SetComm( mpi::Comm comm )
{
    sources_.clear();
    targets_.clear();
    sorted_ = true;
    assembling_ = false;
    localEdgeOffsets_.clear();

    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::CommDup( comm, comm_ );

    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    blocksize_ = numSources_/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank < commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numSources_ - (commSize-1)*blocksize_;
}

inline mpi::Comm 
DistGraph::Comm() const
{ return comm_; }

inline int
DistGraph::Blocksize() const
{ return blocksize_; }

inline int
DistGraph::FirstLocalSource() const
{ return firstLocalSource_; }

inline int
DistGraph::NumLocalSources() const
{ return numLocalSources_; }

inline int
DistGraph::NumLocalEdges() const
{
#ifndef RELEASE
    PushCallStack("DistGraph::NumLocalEdges");
    EnsureConsistentSizes();
    PopCallStack();
#endif
    return sources_.size();
}

inline int
DistGraph::Capacity() const
{
#ifndef RELEASE
    PushCallStack("DistGraph::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
    PopCallStack();
#endif
    return sources_.capacity();
}

inline int
DistGraph::Source( int localEdge ) const
{
#ifndef RELEASE
    PushCallStack("DistGraph::Source");
    if( localEdge < 0 || localEdge >= (int)sources_.size() )
        throw std::logic_error("Edge number out of bounds");
#endif
    EnsureNotAssembling();
#ifndef RELEASE
    PopCallStack();
#endif
    return sources_[localEdge];
}

inline int
DistGraph::Target( int localEdge ) const
{
#ifndef RELEASE
    PushCallStack("DistGraph::Target");
    if( localEdge < 0 || localEdge >= (int)targets_.size() )
        throw std::logic_error("Edge number out of bounds");
#endif
    EnsureNotAssembling();
#ifndef RELEASE
    PopCallStack();
#endif
    return targets_[localEdge];
}

inline int
DistGraph::LocalEdgeOffset( int localSource ) const
{
#ifndef RELEASE
    PushCallStack("DistGraph::LocalEdgeOffset");
    if( localSource < 0 || localSource > numLocalSources_ )
    {
        std::ostringstream msg;
        msg << "Out of bounds localSource: " << localSource << " is not in ["
            << "0," << numLocalSources_ << ")";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    EnsureNotAssembling();
    const int localEdgeOffset = localEdgeOffsets_[localSource];
#ifndef RELEASE
    PopCallStack();
#endif
    return localEdgeOffset;
}

inline int
DistGraph::NumConnections( int localSource ) const
{
#ifndef RELEASE
    PushCallStack("DistGraph::NumConnections");
#endif
    const int numConnections = LocalEdgeOffset(localSource+1) - 
                               LocalEdgeOffset(localSource);
#ifndef RELEASE
    PopCallStack();
#endif
    return numConnections;
}

inline const DistGraph&
DistGraph::operator=( const Graph& graph )
{
#ifndef RELEASE
    PushCallStack("DistGraph::operator=");
#endif
    numSources_ = graph.numSources_; 
    numTargets_ = graph.numTargets_;

    SetComm( mpi::COMM_SELF );

    sources_ = graph.sources_;
    targets_ = graph.targets_;

    sorted_ = graph.sorted_;
    assembling_ = graph.assembling_;
    localEdgeOffsets_ = graph.edgeOffsets_;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

inline const DistGraph& 
DistGraph::operator=( const DistGraph& graph )
{
#ifndef RELEASE
    PushCallStack("DistGraph::operator=");
#endif
    numSources_ = graph.numSources_;
    numTargets_ = graph.numTargets_;

    SetComm( graph.comm_ );

    sources_ = graph.sources_;
    targets_ = graph.targets_;

    sorted_ = graph.sorted_;
    assembling_ = graph.assembling_;
    localEdgeOffsets_ = graph.localEdgeOffsets_;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

inline bool
DistGraph::ComparePairs
( const std::pair<int,int>& a, const std::pair<int,int>& b )
{
    return a.first < b.first || (a.first == b.first && a.second < b.second);
}

inline void
DistGraph::StartAssembly()
{
#ifndef RELEASE
    PushCallStack("DistGraph::StartAssembly");
#endif
    EnsureNotAssembling();
    assembling_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistGraph::StopAssembly()
{
#ifndef RELEASE
    PushCallStack("DistGraph::StopAssembly");
#endif
    if( !assembling_ )
        throw std::logic_error("Cannot stop assembly without starting");
    assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !sorted_ )
    {
        const int numLocalEdges = sources_.size();
        std::vector<std::pair<int,int> > pairs( numLocalEdges );
        for( int e=0; e<numLocalEdges; ++e )
        {
            pairs[e].first = sources_[e];
            pairs[e].second = targets_[e];
        }
        std::sort( pairs.begin(), pairs.end(), ComparePairs );

        // Compress out duplicates
        int lastUnique=0;
        for( int e=1; e<numLocalEdges; ++e )
            if( pairs[e] != pairs[lastUnique] )
                pairs[++lastUnique] = pairs[e];
        const int numUnique = lastUnique+1;

        sources_.resize( numUnique );
        targets_.resize( numUnique );
        for( int e=0; e<numUnique; ++e )
        {
            sources_[e] = pairs[e].first;
            targets_[e] = pairs[e].second;
        }
    }

    ComputeLocalEdgeOffsets();
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistGraph::ComputeLocalEdgeOffsets()
{
#ifndef RELEASE
    PushCallStack("DistGraph::ComputeLocalEdgeOffsets");
#endif
    // Compute the local edge offsets
    int sourceOffset = 0;
    int prevSource = firstLocalSource_-1;
    localEdgeOffsets_.resize( numLocalSources_+1 );
    const int numLocalEdges = NumLocalEdges();
    for( int localEdge=0; localEdge<numLocalEdges; ++localEdge )
    {
        const int source = Source( localEdge );
#ifndef RELEASE
        if( source < prevSource )
            throw std::runtime_error("sources were not properly sorted");
#endif
        while( source != prevSource )
        {
            localEdgeOffsets_[sourceOffset++] = localEdge;
            ++prevSource;
        }
    }
    localEdgeOffsets_[numLocalSources_] = numLocalEdges;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistGraph::Reserve( int numLocalEdges )
{ 
    sources_.reserve( numLocalEdges );
    targets_.reserve( numLocalEdges );
}

inline void
DistGraph::Insert( int source, int target )
{
#ifndef RELEASE
    PushCallStack("DistGraph::Insert");
    EnsureConsistentSizes();
    const int capacity = Capacity();
    const int numLocalEdges = NumLocalEdges();
    if( source < firstLocalSource_ || 
        source >= firstLocalSource_+numLocalSources_ )
    {
        std::ostringstream msg;
        msg << "Source was out of bounds: " << source << " is not in ["
            << firstLocalSource_ << "," << firstLocalSource_+numLocalSources_
            << ")";
        throw std::logic_error( msg.str().c_str() );
    }
    if( numLocalEdges == capacity )
        std::cerr << "WARNING: Pushing back without first reserving space" 
                  << std::endl;
#endif
    if( !assembling_ )
        throw std::logic_error("Must start assembly before pushing back");
    if( sorted_ && sources_.size() != 0 )
    {
        if( source < sources_.back() )
            sorted_ = false;
        if( source == sources_.back() && target < targets_.back() )
            sorted_ = false;
    }
    sources_.push_back( source );
    targets_.push_back( target );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
DistGraph::Empty()
{
    numSources_ = 0;
    numTargets_ = 0;
    sources_.clear();
    targets_.clear();
    blocksize_ = 0;
    firstLocalSource_ = 0;
    numLocalSources_ = 0;
    sorted_ = true;
    assembling_ = false;
    localEdgeOffsets_.clear();
}

inline void
DistGraph::ResizeTo( int numVertices )
{ ResizeTo( numVertices, numVertices ); }

inline void
DistGraph::ResizeTo( int numSources, int numTargets )
{
    const int commRank = mpi::CommRank( comm_ );
    const int commSize = mpi::CommSize( comm_ );
    numSources_ = numSources;
    numTargets_ = numTargets;
    blocksize_ = numSources/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank < commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numSources - (commSize-1)*blocksize_;
    sources_.clear();
    targets_.clear();
    sorted_ = true;
    assembling_ = false;
    localEdgeOffsets_.clear();
}

inline void
DistGraph::EnsureNotAssembling() const
{
    if( assembling_ )
        throw std::logic_error("Should have finished assembling first");
}

inline void
DistGraph::EnsureConsistentSizes() const
{ 
    if( sources_.size() != targets_.size() )
        throw std::logic_error("Inconsistent graph sizes");
}

inline void
DistGraph::EnsureConsistentCapacities() const
{ 
    if( sources_.capacity() != targets_.capacity() )
        throw std::logic_error("Inconsistent graph capacities");
}

} // namespace cliq

#endif // CLIQUE_DIST_GRAPH_MAIN_HPP
