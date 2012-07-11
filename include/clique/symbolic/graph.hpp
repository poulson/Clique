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
#ifndef CLIQUE_GRAPH_HPP
#define CLIQUE_GRAPH_HPP 1

namespace cliq {

class Graph
{
public:
    Graph();
    Graph( int numVertices );
    Graph( int numSources, int numTargets );
    ~Graph();

    int NumSources() const;
    int NumTargets() const;

    int NumEdges() const;
    int Capacity() const;

    int Source( int edge ) const;
    int Target( int edge ) const;

    int EdgeOffset( int source ) const;
    int NumConnections( int source ) const;

    void Reserve( int numEdges );
    void PushBack( int source, int target );

    void Empty();
    void ResizeTo( int numVertices );
    void ResizeTo( int numSources, int numTargets );

private:
    int numSources_, numTargets_;
    std::vector<int> sources_, targets_;

    // Helpers for local indexing
    mutable bool haveEdgeOffsets_;
    mutable std::vector<int> edgeOffsets_;
    void UpdateEdgeOffsets() const;

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    template<typename F> friend class SparseMatrix;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline 
Graph::Graph()
: numSources_(0), numTargets_(0)
{ haveEdgeOffsets_ = false; }

inline 
Graph::Graph( int numVertices )
: numSources_(numVertices), numTargets_(numVertices)
{ haveEdgeOffsets_ = false; }

inline 
Graph::Graph( int numSources, int numTargets )
: numSources_(numSources), numTargets_(numTargets)
{ haveEdgeOffsets_ = false; }

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
    PushCallStack("Graph::NumEdges");
    EnsureConsistentSizes();
    PopCallStack();
#endif
    return sources_.size();
}

inline int
Graph::Capacity() const
{
#ifndef RELEASE
    PushCallStack("Graph::Capacity");
    EnsureConsistentSizes();
    EnsureConsistentCapacities();
    PopCallStack();
#endif
    return sources_.capacity();
}

inline int
Graph::Source( int edge ) const
{
#ifndef RELEASE
    PushCallStack("Graph::Source");
    if( edge < 0 || edge >= sources_.size() )
        throw std::logic_error("Edge number out of bounds");
    PopCallStack();
#endif
    return sources_[edge];
}

inline int
Graph::Target( int edge ) const
{
#ifndef RELEASE
    PushCallStack("Graph::Target");
    if( edge < 0 || edge >= targets_.size() )
        throw std::logic_error("Edge number out of bounds");
    PopCallStack();
#endif
    return targets_[edge];
}

inline int
Graph::EdgeOffset( int source ) const
{
#ifndef RELEASE
    PushCallStack("Graph::EdgeOffset");
    if( source < 0 || source > numSources_ )
        throw std::logic_error("Out of bounds source index");
#endif
    UpdateEdgeOffsets();
    const int edgeOffset = edgeOffsets_[source];
#ifndef RELEASE
    PopCallStack();
#endif
    return edgeOffset;
}

inline int
Graph::NumConnections( int source ) const
{
#ifndef RELEASE
    PushCallStack("Graph::NumConnections");
#endif
    const int numConnections = EdgeOffset(source+1) - EdgeOffset(source);
#ifndef RELEASE
    PopCallStack();
#endif
    return numConnections;
}

inline void
Graph::UpdateEdgeOffsets() const
{
#ifndef RELEASE
    PushCallStack("Graph::UpdateEdgeOffsets");
#endif
    if( !haveEdgeOffsets_ )
    {
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
        haveEdgeOffsets_ = true;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Graph::Reserve( int numEdges )
{ 
    sources_.reserve( numEdges );
    targets_.reserve( numEdges );
}

inline void
Graph::PushBack( int source, int target )
{
#ifndef RELEASE
    PushCallStack("Graph::PushBack");
    EnsureConsistentSizes();
    if( sources_.size() != 0 && source < sources_.back() )
        throw std::logic_error("Incorrectly ordered sources");
    if( targets_.size() != 0 && 
        source == sources_.back() && target < targets_.back() )
        throw std::logic_error("Incorrectly ordered targets");
    const int capacity = Capacity();
    const int numEdges = NumEdges();
    if( numEdges == capacity )
        std::cerr << "WARNING: Pushing back without first reserving space" 
                  << std::endl;
#endif
    sources_.push_back( source );
    targets_.push_back( target );
    haveEdgeOffsets_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Graph::Empty()
{
    numSources_ = 0;
    numTargets_ = 0;
    sources_.clear();
    targets_.clear();
    haveEdgeOffsets_ = false;
    edgeOffsets_.clear();
}

inline void
Graph::ResizeTo( int numVertices )
{ ResizeTo( numVertices, numVertices ); }

inline void
Graph::ResizeTo( int numSources, int numTargets )
{
    numSources_ = numSources;
    numTargets_ = numTargets;
    sources_.clear();
    targets_.clear();
    haveEdgeOffsets_ = false;
    edgeOffsets_.clear();
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

#endif // CLIQUE_GRAPH_HPP
