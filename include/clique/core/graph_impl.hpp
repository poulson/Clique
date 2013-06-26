/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

inline 
Graph::Graph()
: numSources_(0), numTargets_(0), assembling_(false), sorted_(true)
{ }

inline 
Graph::Graph( int numVertices )
: numSources_(numVertices), numTargets_(numVertices), 
  assembling_(false), sorted_(true)
{ }

inline 
Graph::Graph( int numSources, int numTargets )
: numSources_(numSources), numTargets_(numTargets),
  assembling_(false), sorted_(true)
{ }

inline
Graph::Graph( const Graph& graph )
{
#ifndef RELEASE
    CallStackEntry entry("Graph::Graph");
#endif
    if( &graph != this )
        *this = graph;
    else
        throw std::logic_error("Tried to construct a graph with itself");
}
    
inline
Graph::Graph( const DistGraph& graph )
{
#ifndef RELEASE
    CallStackEntry entry("Graph::Graph");
#endif
    *this = graph;
}

inline 
Graph::~Graph()
{ }

inline int 
Graph::NumSources() const
{ return numSources_; }

inline int 
Graph::NumTargets() const
{ return numTargets_; }

inline int
Graph::NumEdges() const
{
#ifndef RELEASE
    CallStackEntry entry("Graph::NumEdges");
    EnsureConsistentSizes();
#endif
    return sources_.size();
}

inline int
Graph::Capacity() const
{
#ifndef RELEASE
    CallStackEntry entry("Graph::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
#endif
    return sources_.capacity();
}

inline int
Graph::Source( int edge ) const
{
#ifndef RELEASE
    CallStackEntry entry("Graph::Source");
    if( edge < 0 || edge >= (int)sources_.size() )
        throw std::logic_error("Edge number out of bounds");
#endif
    EnsureNotAssembling();
    return sources_[edge];
}

inline int
Graph::Target( int edge ) const
{
#ifndef RELEASE
    CallStackEntry entry("Graph::Target");
    if( edge < 0 || edge >= (int)targets_.size() )
        throw std::logic_error("Edge number out of bounds");
#endif
    EnsureNotAssembling();
    return targets_[edge];
}

inline int
Graph::EdgeOffset( int source ) const
{
#ifndef RELEASE
    CallStackEntry entry("Graph::EdgeOffset");
    if( source < 0 )
        throw std::logic_error("Negative source index");
    if( source > numSources_ )
    {
        std::ostringstream msg;
        msg << "Source index was too large: " << source << " is not in "
            << "[0," << numSources_ << "]";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    EnsureNotAssembling();
    return edgeOffsets_[source];
}

inline int
Graph::NumConnections( int source ) const
{
#ifndef RELEASE
    CallStackEntry entry("Graph::NumConnections");
#endif
    return EdgeOffset(source+1) - EdgeOffset(source);
}

inline int*
Graph::SourceBuffer()
{ return &sources_[0]; }

inline int*
Graph::TargetBuffer()
{ return &targets_[0]; }

inline const int*
Graph::LockedSourceBuffer() const
{ return &sources_[0]; }

inline const int*
Graph::LockedTargetBuffer() const
{ return &targets_[0]; }

inline const Graph&
Graph::operator=( const Graph& graph )
{
#ifndef RELEASE
    CallStackEntry entry("Graph::operator=");
#endif
    numSources_ = graph.numSources_;
    numTargets_ = graph.numTargets_;
    sources_ = graph.sources_; 
    targets_ = graph.targets_;

    sorted_ = graph.sorted_;
    assembling_ = graph.assembling_;
    edgeOffsets_ = graph.edgeOffsets_;
    return *this;
}

inline const Graph&
Graph::operator=( const DistGraph& graph )
{
#ifndef RELEASE
    CallStackEntry entry("Graph::operator=");
#endif
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    if( commSize != 1 )
        throw std::logic_error
        ("Cannot yet construct sequential graph from distributed graph");

    numSources_ = graph.numSources_;
    numTargets_ = graph.numTargets_;
    sources_ = graph.sources_; 
    targets_ = graph.targets_;

    sorted_ = graph.sorted_;
    assembling_ = graph.assembling_;
    edgeOffsets_ = graph.localEdgeOffsets_;
    return *this;
}

inline bool
Graph::ComparePairs
( const std::pair<int,int>& a, const std::pair<int,int>& b )
{ return a.first < b.first || (a.first  == b.first && a.second < b.second); }

inline void
Graph::StartAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("Graph::StartAssembly");
#endif
    EnsureNotAssembling();
    assembling_ = true;
}

inline void
Graph::StopAssembly()
{
#ifndef RELEASE
    CallStackEntry entry("Graph::StopAssembly");
#endif
    if( !assembling_ )
        throw std::logic_error("Cannot stop assembly without starting");
    assembling_ = false;

    // Ensure that the connection pairs are sorted
    if( !sorted_ )
    {
        const int numEdges = sources_.size();
        std::vector<std::pair<int,int> > pairs( numEdges );
        for( int e=0; e<numEdges; ++e )
        {
            pairs[e].first = sources_[e];
            pairs[e].second = targets_[e];
        }
        std::sort( pairs.begin(), pairs.end(), ComparePairs );

        // Compress out duplicates
        int lastUnique=0;
        for( int e=1; e<numEdges; ++e )
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

    ComputeEdgeOffsets();
}

inline void
Graph::ComputeEdgeOffsets()
{
#ifndef RELEASE
    CallStackEntry entry("Graph::ComputeEdgeOffsets");
#endif
    // Compute the edge offsets
    int sourceOffset = 0;
    int prevSource = -1;
    edgeOffsets_.resize( numSources_+1 );
    const int numEdges = NumEdges();
    for( int edge=0; edge<numEdges; ++edge )
    {
        const int source = Source( edge );
#ifndef RELEASE
        if( source < prevSource )
            throw std::runtime_error("sources were not properly sorted");
#endif
        while( source != prevSource )
        {
            edgeOffsets_[sourceOffset++] = edge;
            ++prevSource;
        }
    }
    edgeOffsets_[numSources_] = numEdges;
}

inline void
Graph::Reserve( int numEdges )
{ 
    sources_.reserve( numEdges );
    targets_.reserve( numEdges );
}

inline void
Graph::Insert( int source, int target )
{
#ifndef RELEASE
    CallStackEntry entry("Graph::Insert");
    EnsureConsistentSizes();
    const int capacity = Capacity();
    const int numEdges = NumEdges();
    if( source < 0 || source >= numSources_ )
    {
        std::ostringstream msg;
        msg << "Source was out of bounds: " << source << " is not in [0,"
            << numSources_ << ")";
        throw std::logic_error( msg.str().c_str() );
    }
    if( numEdges == capacity )
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
}

inline void
Graph::Empty()
{
    numSources_ = 0;
    numTargets_ = 0;
    std::vector<int>().swap( sources_ );
    std::vector<int>().swap( targets_ );
    sorted_ = true;
    assembling_ = false;
    std::vector<int>().swap( edgeOffsets_ );
}

inline void
Graph::ResizeTo( int numVertices )
{ ResizeTo( numVertices, numVertices ); }

inline void
Graph::ResizeTo( int numSources, int numTargets )
{
    numSources_ = numSources;
    numTargets_ = numTargets;
    std::vector<int>().swap( sources_ );
    std::vector<int>().swap( targets_ );
    sorted_ = true;
    assembling_ = false;
    std::vector<int>().swap( edgeOffsets_ );
}

inline void
Graph::EnsureNotAssembling() const
{
    if( assembling_ )
        throw std::logic_error("Should have finished assembling first");
}

inline void
Graph::EnsureConsistentSizes() const
{ 
    if( sources_.size() != targets_.size() )
        throw std::logic_error("Inconsistent graph sizes");
}

inline void
Graph::EnsureConsistentCapacities() const
{ 
    if( sources_.capacity() != targets_.capacity() )
        throw std::logic_error("Inconsistent graph capacities");
}

} // namespace cliq
