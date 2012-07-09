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
#ifndef CLIQUE_DIST_GRAPH_HPP
#define CLIQUE_DIST_GRAPH_HPP 1

namespace cliq {

// Use a simple 1d distribution where each process owns a fixed number of 
// sources:
//     if last process,  numSources - (commSize-1)*ceil(numSources/commSize)
//     otherwise,        ceil(numSources/commSize)
class DistGraph
{
public:
    DistGraph( int numVertices, mpi::Comm comm );
    DistGraph( int numSources, int numTargets, mpi::Comm comm );

    int NumSources() const;
    int NumTargets() const;

    mpi::Comm Comm() const;
    int Blocksize() const;
    int FirstLocalSource() const;
    int NumLocalSources() const;

    int NumLocalEdges() const;
    int Capacity() const;

    void Reserve( int numLocalEdges );
    void PushBack( int i, int j );

    void Empty();
    void ResizeTo( int numVertices );
    void ResizeTo( int numSources, int numTargets );

private:
    int numSources_, numTargets_;
    mpi::Comm comm_;

    int blocksize_;
    int firstLocalSource_, numLocalSources_;

    std::vector<int> sources_, targets_;

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    template<typename F> friend class DistSparseMatrix;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline DistGraph::DistGraph
( int numVertices, mpi::Comm comm )
: numSources_(numVertices), numTargets_(numVertices), comm_(comm)
{
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    blocksize_ = (numVertices+commSize-1)/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank != commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numVertices - (commSize-1)*blocksize_;
}

inline DistGraph::DistGraph
( int numSources, int numTargets, mpi::Comm comm )
: numSources_(numSources), numTargets_(numTargets), comm_(comm)
{
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    blocksize_ = (numSources+commSize-1)/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank != commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numSources - (commSize-1)*blocksize_;
}

inline int 
DistGraph::NumSources() const
{ return numSources_; }

inline int 
DistGraph::NumTargets() const
{ return numTargets_; }

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

inline void
DistGraph::Reserve( int numLocalEdges )
{ 
    sources_.reserve( numLocalEdges );
    targets_.reserve( numLocalEdges );
}

inline void
DistGraph::PushBack( int source, int target )
{
#ifndef RELEASE
    PushCallStack("DistGraph::PushBack");
    EnsureConsistentSizes();
    if( sources_.size() != 0 && source < sources_.back() )
        throw std::logic_error("Incorrectly ordered sources");
    if( targets_.size() != 0 && target < targets_.back() )
        throw std::logic_error("Incorrectly ordered targets");
    const int capacity = Capacity();
    const int numLocalEdges = NumLocalEdges();
    if( numLocalEdges == capacity )
        std::cerr << "WARNING: Pushing back without first reserving space" 
                  << std::endl;
#endif
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
    blocksize_ = (numSources+commSize-1)/commSize;
    firstLocalSource_ = commRank*blocksize_;
    if( commRank != commSize-1 )
        numLocalSources_ = blocksize_;
    else
        numLocalSources_ = numSources - (commSize-1)*blocksize_;
    sources_.clear();
    targets_.clear();
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

#endif // CLIQUE_DIST_GRAPH_HPP
