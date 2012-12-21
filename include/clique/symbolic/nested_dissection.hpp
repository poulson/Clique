/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_PARMETIS
# include "parmetis.h"

extern "C" {

void CliqBisect
( idx_t* nvtxs, idx_t* xAdj, idx_t* adjacency, idx_t* nseps, real_t* imbalance,
  idx_t* perm, idx_t* sizes );

void CliqParallelBisect
( idx_t* vtxDist, idx_t* xAdj, idx_t* adjacency, 
  idx_t* nparseps, idx_t* nseqseps, real_t* imbalance, idx_t* options, 
  idx_t* perm, idx_t* sizes, MPI_Comm* comm );

} // extern "C"
#endif 

namespace cliq {

#ifdef HAVE_PARMETIS
void NestedDissection
( const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
        bool sequential=true,
        int numDistSeps=1, 
        int numSeqSeps=1, 
        int cutoff=128, 
        bool storeFactRecvIndices=false );

int Bisect
( const Graph& graph, 
        Graph& leftChild, 
        Graph& rightChild, 
        std::vector<int>& perm, 
        int numSeps=5 );

// NOTE: for two or more processes
int Bisect
( const DistGraph& graph, 
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft,
        bool sequential=true,
        int numDistSeps=1, 
        int numSeqSeps=1 );
#endif // HAVE_PARMETIS

int DistributedDepth( mpi::Comm comm );
void EnsurePermutation( const std::vector<int>& map );
void EnsurePermutation( const DistMap& map );
void ReverseOrder( DistSeparatorTree& sepTree, DistSymmElimTree& eTree );

void BuildChildrenFromPerm
( const Graph& graph, const std::vector<int>& perm, 
  int leftChildSize, Graph& leftChild,
  int rightChildSize, Graph& rightChild );
void BuildChildFromPerm
( const DistGraph& graph, const DistMap& perm,
  int leftChildSize, int rightChildSize,
  bool& onLeft, DistGraph& child );

void BuildMap
( const DistGraph& graph, 
  const DistSeparatorTree& sepTree, 
        DistMap& map );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline void
DistributedDepthRecursion
( unsigned commRank, unsigned commSize, unsigned& distDepth )
{
    if( commSize == 1 )
        return;

    ++distDepth;
    const unsigned smallTeamSize = commSize/2;
    const unsigned largeTeamSize = commSize - smallTeamSize;
    if( commRank < smallTeamSize )
        DistributedDepthRecursion( commRank, smallTeamSize, distDepth );
    else
        DistributedDepthRecursion
        ( commRank-smallTeamSize, largeTeamSize, distDepth );
}

inline int
DistributedDepth( mpi::Comm comm )
{
    unsigned commRank = mpi::CommRank( comm );
    unsigned commSize = mpi::CommSize( comm );
    unsigned distDepth = 0;
    DistributedDepthRecursion( commRank, commSize, distDepth );
    return distDepth;
}

#ifdef HAVE_PARMETIS
inline void
NestedDissectionRecursion
( const Graph& graph, 
  const std::vector<int>& perm,
        DistSeparatorTree& sepTree, 
        DistSymmElimTree& eTree,
        int parent, 
        int offset, 
        int numSeps=5,
        int cutoff=128 )
{
#ifndef RELEASE
    PushCallStack("NestedDissectionRecursion");
#endif
    if( graph.NumSources() <= cutoff )
    {
        // Fill in this node of the local separator tree
        const int numSources = graph.NumSources();
        sepTree.localSepsAndLeaves.push_back( new LocalSepOrLeaf );
        LocalSepOrLeaf& leaf = *sepTree.localSepsAndLeaves.back();
        leaf.parent = parent;
        leaf.offset = offset;
        leaf.indices = perm;

        // Fill in this node of the local elimination tree
        eTree.localNodes.push_back( new LocalSymmNode );
        LocalSymmNode& node = *eTree.localNodes.back();
        node.size = numSources;
        node.offset = offset;
        node.parent = parent;
        node.children.clear();
        std::set<int> connectedAncestors;
        for( int s=0; s<node.size; ++s )
        {
            const int numConnections = graph.NumConnections( s );
            const int edgeOffset = graph.EdgeOffset( s );
            for( int t=0; t<numConnections; ++t )
            {
                const int target = graph.Target( edgeOffset+t );
                if( target >= numSources )
                    connectedAncestors.insert( offset+target );
            }
        }
        node.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          node.lowerStruct.begin() );
    }
    else
    {
        // Partition the graph and construct the inverse map
        Graph leftChild, rightChild;
        std::vector<int> map;
        const int sepSize = 
            Bisect( graph, leftChild, rightChild, map, numSeps );
        const int numSources = graph.NumSources();
        std::vector<int> inverseMap( numSources );
        for( int s=0; s<numSources; ++s )
            inverseMap[map[s]] = s;

        // Mostly compute this node of the local separator tree
        // (we will finish computing the separator indices soon)
        sepTree.localSepsAndLeaves.push_back( new LocalSepOrLeaf );
        LocalSepOrLeaf& sep = *sepTree.localSepsAndLeaves.back();
        sep.parent = parent;
        sep.offset = offset + (numSources-sepSize);
        sep.indices.resize( sepSize );
        for( int s=0; s<sepSize; ++s )
        {
            const int mappedSource = s + (numSources-sepSize);
            sep.indices[s] = inverseMap[mappedSource];
        }
    
        // Fill in this node in the local elimination tree
        eTree.localNodes.push_back( new LocalSymmNode );
        LocalSymmNode& node = *eTree.localNodes.back();
        node.size = sepSize;
        node.offset = sep.offset;
        node.parent = parent;
        node.children.resize( 2 );
        std::set<int> connectedAncestors;
        for( int s=0; s<sepSize; ++s )
        {
            const int source = sep.indices[s];
            const int numConnections = graph.NumConnections( source );
            const int edgeOffset = graph.EdgeOffset( source );
            for( int t=0; t<numConnections; ++t )
            {
                const int target = graph.Target( edgeOffset+t );
                if( target >= numSources )
                    connectedAncestors.insert( offset+target );
            }
        }
        node.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          node.lowerStruct.begin() );

        // Finish computing the separator indices
        for( int s=0; s<sepSize; ++s )
            sep.indices[s] = perm[sep.indices[s]];

        // Construct the inverse maps from the child indices to the original
        // degrees of freedom
        const int leftChildSize = leftChild.NumSources();
        std::vector<int> leftPerm( leftChildSize );
        for( int s=0; s<leftChildSize; ++s )
            leftPerm[s] = perm[inverseMap[s]];
        const int rightChildSize = rightChild.NumSources();
        std::vector<int> rightPerm( rightChildSize );
        for( int s=0; s<rightChildSize; ++s )
            rightPerm[s] = perm[inverseMap[s+leftChildSize]];

        // Update right then left so that, once we later reverse the order 
        // of the nodes, the left node will be ordered first
        const int parent = eTree.localNodes.size()-1;
        node.children[1] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( rightChild, rightPerm, sepTree, eTree, parent, offset+leftChildSize,
          numSeps, cutoff );
        node.children[0] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( leftChild, leftPerm, sepTree, eTree, parent, offset, 
          numSeps, cutoff );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
NestedDissectionRecursion
( const DistGraph& graph, 
  const DistMap& perm,
        DistSeparatorTree& sepTree, 
        DistSymmElimTree& eTree,
        int depth, 
        int offset, 
        bool onLeft,
        bool sequential=true,
        int numDistSeps=1, 
        int numSeqSeps=1,
        int cutoff=128 )
{
#ifndef RELEASE
    PushCallStack("NestedDissectionRecursion");
#endif
    const int distDepth = sepTree.distSeps.size();
    mpi::Comm comm = graph.Comm();
    if( distDepth - depth > 0 )
    {
        // Partition the graph and construct the inverse map
        DistGraph child;
        bool childIsOnLeft;
        DistMap map;
        const int sepSize = 
            Bisect
            ( graph, child, map, childIsOnLeft, 
              sequential, numDistSeps, numSeqSeps );
        const int numSources = graph.NumSources();
        const int childSize = child.NumSources();
        const int leftChildSize = 
            ( childIsOnLeft ? childSize : numSources-sepSize-childSize );

        DistMap inverseMap;
        map.FormInverse( inverseMap );

        // Mostly fill this node of the DistSeparatorTree
        // (we will finish computing the separator indices at the end)
        DistSeparator& sep = sepTree.distSeps[distDepth-1-depth];
        mpi::CommDup( comm, sep.comm );
        sep.offset = offset + (numSources-sepSize);
        sep.indices.resize( sepSize );
        for( int s=0; s<sepSize; ++s )
            sep.indices[s] = s + (numSources-sepSize);
        inverseMap.Translate( sep.indices );

        // Fill in this node of the DistSymmElimTree
        DistSymmNode& node = eTree.distNodes[distDepth-depth];
        node.size = sepSize;
        node.offset = sep.offset;
        node.onLeft = onLeft;
        mpi::CommDup( comm, node.comm );
        const int numLocalSources = graph.NumLocalSources();
        const int firstLocalSource = graph.FirstLocalSource();
        std::set<int> localConnectedAncestors;
        for( int s=0; s<sepSize; ++s )
        {
            const int source = sep.indices[s];
            if( source >= firstLocalSource && 
                source < firstLocalSource+numLocalSources )
            {
                const int localSource = source - firstLocalSource;
                const int numConnections = graph.NumConnections( localSource );
                const int localOffset = graph.LocalEdgeOffset( localSource );
                for( int t=0; t<numConnections; ++t )
                {
                    const int target = graph.Target( localOffset+t );
                    if( target >= numSources )
                        localConnectedAncestors.insert( offset+target );
                }
            }
        }
        const int numLocalConnected = localConnectedAncestors.size();
        const int commSize = mpi::CommSize( comm );
        std::vector<int> localConnectedSizes( commSize );
        mpi::AllGather
        ( &numLocalConnected, 1, &localConnectedSizes[0], 1, comm );
        std::vector<int> localConnectedVector( numLocalConnected );
        std::copy
        ( localConnectedAncestors.begin(), localConnectedAncestors.end(), 
          localConnectedVector.begin() );
        int sumOfLocalConnectedSizes=0;
        std::vector<int> localConnectedOffsets( commSize );
        for( int q=0; q<commSize; ++q )
        {
            localConnectedOffsets[q] = sumOfLocalConnectedSizes;
            sumOfLocalConnectedSizes += localConnectedSizes[q];
        }
        std::vector<int> localConnections( sumOfLocalConnectedSizes );
        mpi::AllGather
        ( &localConnectedVector[0], numLocalConnected,
          &localConnections[0], 
          &localConnectedSizes[0], &localConnectedOffsets[0], comm );
        std::set<int> connectedAncestors
        ( localConnections.begin(), localConnections.end() );
        node.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          node.lowerStruct.begin() );

        // Finish computing the separator indices
        perm.Translate( sep.indices );

        // Construct map from child indices to the original ordering
        DistMap newPerm( child.NumSources(), child.Comm() );
        const int localChildSize = child.NumLocalSources();
        const int firstLocalChildSource = child.FirstLocalSource();
        if( childIsOnLeft )
            for( int s=0; s<localChildSize; ++s )
                newPerm.SetLocal( s, s+firstLocalChildSource );
        else
            for( int s=0; s<localChildSize; ++s )
                newPerm.SetLocal( s, s+firstLocalChildSource+leftChildSize );
        inverseMap.Extend( newPerm );
        perm.Extend( newPerm );

        // Recurse
        const int newOffset = ( childIsOnLeft ? offset : offset+leftChildSize );
        NestedDissectionRecursion
        ( child, newPerm, sepTree, eTree, depth+1, newOffset, 
          childIsOnLeft, sequential, numDistSeps, numSeqSeps, cutoff );
    }
    else if( graph.NumSources() <= cutoff )
    {
        // Convert to a sequential graph
        const int numSources = graph.NumSources();
        Graph seqGraph( graph );

        // Fill in this node of the local separator tree
        sepTree.localSepsAndLeaves.push_back( new LocalSepOrLeaf );
        LocalSepOrLeaf& leaf = *sepTree.localSepsAndLeaves.back();
        leaf.parent = -1;
        leaf.offset = offset;
        leaf.indices = perm.LocalMap();

        // Fill in this node of the local and distributed parts of the 
        // elimination tree
        eTree.localNodes.push_back( new LocalSymmNode );
        LocalSymmNode& localNode = *eTree.localNodes.back();
        DistSymmNode& distNode = eTree.distNodes[0];
        mpi::CommDup( comm, distNode.comm );
        distNode.onLeft = onLeft;
        distNode.size = localNode.size = numSources;
        distNode.offset = localNode.offset = offset;
        localNode.parent = -1;
        localNode.children.clear();
        std::set<int> connectedAncestors;
        for( int s=0; s<numSources; ++s )
        {
            const int numConnections = seqGraph.NumConnections( s );
            const int edgeOffset = seqGraph.EdgeOffset( s );
            for( int t=0; t<numConnections; ++t )
            {
                const int target = seqGraph.Target( edgeOffset+t );
                if( target >= numSources )
                    connectedAncestors.insert( offset+target );
            }
        }
        localNode.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          localNode.lowerStruct.begin() );    
        distNode.lowerStruct = localNode.lowerStruct;
    }
    else
    {
        // Convert to a sequential graph
        Graph seqGraph( graph );

        // Partition the graph and construct the inverse map
        Graph leftChild, rightChild;
        std::vector<int> map;
        const int sepSize = 
            Bisect( seqGraph, leftChild, rightChild, map, numSeqSeps );
        const int numSources = graph.NumSources();
        std::vector<int> inverseMap( numSources );
        for( int s=0; s<numSources; ++s )
            inverseMap[map[s]] = s;

        // Mostly compute this node of the local separator tree
        // (we will finish computing the separator indices soon)
        sepTree.localSepsAndLeaves.push_back( new LocalSepOrLeaf );
        LocalSepOrLeaf& sep = *sepTree.localSepsAndLeaves.back();
        sep.parent = -1;
        sep.offset = offset + (numSources-sepSize);
        sep.indices.resize( sepSize );
        for( int s=0; s<sepSize; ++s )
        {
            const int mappedSource = s + (numSources-sepSize);
            sep.indices[s] = inverseMap[mappedSource];
        }
        
        // Fill in this node in both the local and distributed parts of 
        // the elimination tree
        eTree.localNodes.push_back( new LocalSymmNode );
        LocalSymmNode& localNode = *eTree.localNodes.back();
        DistSymmNode& distNode = eTree.distNodes[0];
        mpi::CommDup( comm, distNode.comm );
        distNode.onLeft = onLeft;
        distNode.size = localNode.size = sepSize;
        distNode.offset = localNode.offset = sep.offset;
        localNode.parent = -1;
        localNode.children.resize( 2 );
        std::set<int> connectedAncestors;
        for( int s=0; s<sepSize; ++s )
        {
            const int source = sep.indices[s];
            const int numConnections = seqGraph.NumConnections( source );
            const int edgeOffset = seqGraph.EdgeOffset( source );
            for( int t=0; t<numConnections; ++t )
            {
                const int target = seqGraph.Target( edgeOffset+t );
                if( target >= numSources )
                    connectedAncestors.insert( offset+target );
            }
        }
        localNode.lowerStruct.resize( connectedAncestors.size() );
        std::copy
        ( connectedAncestors.begin(), connectedAncestors.end(), 
          localNode.lowerStruct.begin() );
        distNode.lowerStruct = localNode.lowerStruct;

        // Finish computing the separator indices
        // (This is a faster version of the Translate member function)
        for( int s=0; s<sepSize; ++s )
            sep.indices[s] = perm.GetLocal( sep.indices[s] );

        // Construct the inverse maps from the child indices to the original
        // degrees of freedom
        const int leftChildSize = leftChild.NumSources();
        std::vector<int> leftPerm( leftChildSize );
        for( int s=0; s<leftChildSize; ++s )
            leftPerm[s] = perm.GetLocal( inverseMap[s] );
        const int rightChildSize = rightChild.NumSources();
        std::vector<int> rightPerm( rightChildSize );
        for( int s=0; s<rightChildSize; ++s )
            rightPerm[s] = perm.GetLocal( inverseMap[s+leftChildSize] );

        // Update right then left so that, once we later reverse the order 
        // of the nodes, the left node will be ordered first
        const int parent=0;
        localNode.children[1] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( rightChild, rightPerm, sepTree, eTree, parent, offset+leftChildSize, 
          numSeqSeps, cutoff );
        localNode.children[0] = eTree.localNodes.size();
        NestedDissectionRecursion
        ( leftChild, leftPerm, sepTree, eTree, parent, offset, 
          numSeqSeps, cutoff );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void 
NestedDissection
( const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
        bool sequential,
        int numDistSeps, 
        int numSeqSeps, 
        int cutoff,
        bool storeFactRecvIndices )
{
#ifndef RELEASE
    PushCallStack("NestedDissection");
#endif
    // NOTE: There is a potential memory leak here if these data structures 
    //       are reused. Their destructors should call a member function which
    //       we can simply call here to clear the data
    DistSymmElimTree eTree;
    eTree.localNodes.clear();
    sepTree.localSepsAndLeaves.clear();

    mpi::Comm comm = graph.Comm();
    const int distDepth = DistributedDepth( comm );
    eTree.distNodes.resize( distDepth+1 );
    sepTree.distSeps.resize( distDepth );

    DistMap perm( graph.NumSources(), graph.Comm() );
    const int firstLocalSource = perm.FirstLocalSource();
    const int numLocalSources = perm.NumLocalSources();
    for( int s=0; s<numLocalSources; ++s )
        perm.SetLocal( s, s+firstLocalSource );
    NestedDissectionRecursion
    ( graph, perm, sepTree, eTree, 0, 0, false, sequential, 
      numDistSeps, numSeqSeps, cutoff );

    ReverseOrder( sepTree, eTree );

    // Construct the distributed reordering    
    BuildMap( graph, sepTree, map );
#ifndef RELEASE
    EnsurePermutation( map );
#endif

    // Run the symbolic analysis
    SymmetricAnalysis( eTree, info, storeFactRecvIndices );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline int 
Bisect
( const Graph& graph ,Graph& leftChild, Graph& rightChild,
  std::vector<int>& perm, int numSeps )
{
#ifndef RELEASE
    PushCallStack("Bisect");
#endif
    // METIS assumes that there are no self-connections or connections 
    // outside the sources, so we must manually remove them from our graph
    const int numSources = graph.NumSources();
    const int numEdges = graph.NumEdges();
    int numValidEdges = 0;
    for( int i=0; i<numEdges; ++i )
        if( graph.Source(i) != graph.Target(i) && graph.Target(i) < numSources )
            ++numValidEdges;

    // Fill our connectivity (ignoring self and too-large connections)
    std::vector<idx_t> xAdj( numSources+1 );
    std::vector<idx_t> adjacency( numValidEdges );
    int validCounter=0;
    int sourceOffset=0;
    int prevSource=-1;
    for( int edge=0; edge<numEdges; ++edge )
    {
        const int source = graph.Source( edge );
        const int target = graph.Target( edge );
#ifndef RELEASE
        if( source < prevSource )
            throw std::runtime_error("sources were not properly sorted");
#endif
        while( source != prevSource )
        {
            xAdj[sourceOffset++] = validCounter;
            ++prevSource;
        }
        if( source != target && target < numSources )
        {
            adjacency[validCounter] = target;
            ++validCounter;
        }
    }
#ifndef RELEASE
    if( sourceOffset != numSources )
        throw std::logic_error("Mistake in xAdj computation");
#endif
    xAdj[numSources] = numValidEdges;

    // Create space for the result
    perm.resize( numSources );

    // Use the custom METIS interface
    idx_t nvtxs = numSources;
    idx_t nseps = numSeps;
    real_t imbalance = 1.1;
    std::vector<idx_t> sizes(3);
    CliqBisect
    ( &nvtxs, &xAdj[0], &adjacency[0], &nseps, &imbalance, &perm[0], 
      &sizes[0] );
#ifndef RELEASE
    EnsurePermutation( perm );
#endif
    BuildChildrenFromPerm
    ( graph, perm, sizes[0], leftChild, sizes[1], rightChild );
#ifndef RELEASE
    PopCallStack();
#endif
    return sizes[2];
}

inline int 
Bisect
( const DistGraph& graph, 
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft, 
        bool sequential,
        int numDistSeps, 
        int numSeqSeps )
{
#ifndef RELEASE
    PushCallStack("Bisect");
#endif
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );
    if( commSize == 1 )
        throw std::logic_error
        ("This routine assumes at least two processes are used, "
         "otherwise one child will be lost");

    // ParMETIS assumes that there are no self-connections or connections 
    // outside the sources, so we must manually remove them from our graph
    const int numSources = graph.NumSources();
    const int numLocalEdges = graph.NumLocalEdges();
    int numLocalValidEdges = 0;
    for( int i=0; i<numLocalEdges; ++i )
        if( graph.Source(i) != graph.Target(i) && graph.Target(i) < numSources )
            ++numLocalValidEdges;

    // Fill our local connectivity (ignoring self and too-large connections)
    const int blocksize = graph.Blocksize();
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
        if( source != target && target < numSources )
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

    idx_t nparseps = numDistSeps;
    idx_t nseqseps = numSeqSeps;
    real_t imbalance = 1.1;
    std::vector<idx_t> sizes(3);
    if( sequential )
    {
        // Gather the number of local valid edges on the root process
        std::vector<int> edgeSizes( commSize ), edgeOffsets;
        mpi::AllGather( &numLocalValidEdges, 1, &edgeSizes[0], 1, comm );
        int numEdges;
        if( commRank == 0 )
        {
            edgeOffsets.resize( commSize );
            numEdges=0;
            for( int q=0; q<commSize; ++q )
            {
                edgeOffsets[q] = numEdges;
                numEdges += edgeSizes[q];
            }
        }

        // Gather the edges on the root process (with padding)
        int maxLocalValidEdges=0;
        for( int q=0; q<commSize; ++q )
            maxLocalValidEdges = 
                std::max( maxLocalValidEdges, edgeSizes[q] );
        adjacency.resize( maxLocalValidEdges );
        std::vector<int> globalAdj;
        if( commRank == 0 )
            globalAdj.resize( maxLocalValidEdges*commSize, 0 );
        mpi::Gather
        ( &adjacency[0], maxLocalValidEdges, 
          &globalAdj[0], maxLocalValidEdges, 0, comm );

        if( commRank == 0 )
        {
            // Remove the padding
            for( int q=1; q<commSize; ++q )
            {
                const int edgeOffset = q*maxLocalValidEdges;
                for( int j=0; j<edgeSizes[q]; ++j )
                    globalAdj[edgeOffsets[q]+j] = globalAdj[edgeOffset+j];
            }
            globalAdj.resize( numEdges );
        }

        // Set up the global xAdj vector
        std::vector<int> globalXAdj;
        // Set up the first commSize*blocksize entries
        if( commRank == 0 )
            globalXAdj.resize( numSources+1 );
        mpi::Gather( &xAdj[0], blocksize, &globalXAdj[0], blocksize, 0, comm );
        if( commRank == 0 )
            for( int q=1; q<commSize; ++q )
                for( int j=0; j<blocksize; ++j )
                    globalXAdj[q*blocksize+j] += edgeOffsets[q];
        // Fix the remaining entries
        const int numRemaining = numSources - commSize*blocksize;
        if( commRank == commSize-1 )
            mpi::Send( &xAdj[blocksize], numRemaining, 0, 0, comm );
        if( commRank == 0 )
        {
            mpi::Recv
            ( &globalXAdj[commSize*blocksize], numRemaining, 
              commSize-1, 0, comm );
            for( int j=0; j<numRemaining; ++j )
                globalXAdj[commSize*blocksize+j] += edgeOffsets[commSize-1];
            globalXAdj[numSources] = numEdges;
        }

        std::vector<int> seqPerm;
        if( commRank == 0 )
        {
            // Use the custom METIS interface
            idx_t nvtxs = numSources;
            seqPerm.resize( nvtxs );
            CliqBisect
            ( &nvtxs, &globalXAdj[0], &globalAdj[0], &nseqseps,
              &imbalance, &seqPerm[0], &sizes[0] );
        }

        // Set up space for the distributed permutation
        perm.SetComm( comm );
        perm.ResizeTo( numSources );

        // Distribute the first commSize*blocksize values of the permutation
        mpi::Scatter
        ( &seqPerm[0], blocksize, perm.LocalBuffer(), blocksize, 0, comm );

        // Make sure the last process gets the straggling entries
        if( commRank == 0 )
            mpi::Send
            ( &seqPerm[commSize*blocksize], numRemaining, commSize-1, 0, comm );
        if( commRank == commSize-1 )
            mpi::Recv
            ( perm.LocalBuffer()+blocksize, numLocalSources-blocksize, 0, 
              0, comm );

        // Broadcast the sizes information from the root
        mpi::Broadcast( (byte*)&sizes[0], 3*sizeof(idx_t), 0, comm );
    }
    else
    {
        // Describe the source distribution
        std::vector<idx_t> vtxDist( commSize+1 );
        for( int i=0; i<commSize; ++i )
            vtxDist[i] = i*blocksize;
        vtxDist[commSize] = graph.NumSources();

        // Create space for the result
        perm.SetComm( comm );
        perm.ResizeTo( numSources );

        // Use the custom ParMETIS interface
        CliqParallelBisect
        ( &vtxDist[0], &xAdj[0], &adjacency[0], &nparseps, &nseqseps, 
          &imbalance, NULL, perm.LocalBuffer(), &sizes[0], &comm );
    }
#ifndef RELEASE
    EnsurePermutation( perm );
#endif
    BuildChildFromPerm( graph, perm, sizes[0], sizes[1], onLeft, child );
#ifndef RELEASE
    PopCallStack();
#endif
    return sizes[2];
}
#endif // HAVE_PARMETIS

inline void
EnsurePermutation( const std::vector<int>& map )
{
#ifndef RELEASE
    PushCallStack("EnsurePermutation");
#endif
    const int numSources = map.size();
    std::vector<int> timesMapped( numSources, 0 );
    for( int i=0; i<numSources; ++i )
        ++timesMapped[map[i]];
    for( int i=0; i<numSources; ++i )
    {
        if( timesMapped[i] != 1 )
        {
            std::ostringstream msg;
            msg << timesMapped[i] << " vertices were relabeled as "
                << i << " in sequential map";
            throw std::logic_error( msg.str().c_str() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
EnsurePermutation( const DistMap& map )
{
#ifndef RELEASE
    PushCallStack("EnsurePermutation");
#endif
    mpi::Comm comm = map.Comm();
    const int commRank = mpi::CommRank( comm );
    const int numSources = map.NumSources();
    const int numLocalSources = map.NumLocalSources();
    std::vector<int> timesMapped( numSources, 0 );
    for( int iLocal=0; iLocal<numLocalSources; ++iLocal )
        ++timesMapped[map.GetLocal(iLocal)];
    mpi::Reduce( &timesMapped[0], numSources, MPI_SUM, 0, comm );
    if( commRank == 0 )
    {
        for( int i=0; i<numSources; ++i )
        {
            if( timesMapped[i] != 1 )
            {
                std::ostringstream msg;
                msg << timesMapped[i] << " vertices were relabeled as "
                    << i << " in parallel map";
                throw std::logic_error( msg.str().c_str() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
ReverseOrder( DistSeparatorTree& sepTree, DistSymmElimTree& eTree )
{
    // Reverse the order of the pointers and indices in the elimination and 
    // separator trees (so that the leaves come first)
    const int numLocalNodes = eTree.localNodes.size();
    const int lastIndex = numLocalNodes-1;
    if( numLocalNodes != 1 )
    {
        // Switch the pointers for the root and last nodes
        LocalSymmNode* rootNode = eTree.localNodes[0];
        LocalSymmNode* lastNode = eTree.localNodes.back();
        eTree.localNodes[0] = lastNode;
        eTree.localNodes.back() = rootNode;
        LocalSepOrLeaf* rootSep = sepTree.localSepsAndLeaves[0];
        LocalSepOrLeaf* lastLeaf = sepTree.localSepsAndLeaves.back();
        sepTree.localSepsAndLeaves[0] = lastLeaf;
        sepTree.localSepsAndLeaves.back() = rootSep;

        // Update their parent indices to what their final values will be
        // (The root node's parent index does not need to be changed.)
        lastNode->parent = lastIndex - lastNode->parent;
        lastLeaf->parent = lastIndex - lastLeaf->parent;
        // Update their children's indices
        // (The last node will not have children)
        const int numRootChildren = rootNode->children.size();
        for( int c=0; c<numRootChildren; ++c )
            rootNode->children[c] = lastIndex - rootNode->children[c];
    }
    // Switch the middle nodes (we will miss the middle node if an odd number)
    for( int s=1; s<numLocalNodes/2; ++s )
    {
        const int t = lastIndex - s;
        // Switch the pointers for the last and right nodes
        LocalSymmNode* leftNode = eTree.localNodes[s];
        LocalSymmNode* rightNode = eTree.localNodes[t];
        LocalSepOrLeaf* leftSepOrLeaf = sepTree.localSepsAndLeaves[s];
        LocalSepOrLeaf* rightSepOrLeaf = sepTree.localSepsAndLeaves[t];
        eTree.localNodes[s] = rightNode;
        eTree.localNodes[t] = leftNode;
        sepTree.localSepsAndLeaves[s] = rightSepOrLeaf;
        sepTree.localSepsAndLeaves[t] = leftSepOrLeaf;
        // Update their parent indices to what their final values will be
        leftNode->parent = lastIndex - leftNode->parent;
        rightNode->parent = lastIndex - rightNode->parent;
        leftSepOrLeaf->parent = lastIndex - leftSepOrLeaf->parent;
        rightSepOrLeaf->parent = lastIndex - rightSepOrLeaf->parent;
        // Update their children's indices
        const int numLeftChildren = leftNode->children.size();
        for( int c=0; c<numLeftChildren; ++c )
            leftNode->children[c] = lastIndex - leftNode->children[c];
        const int numRightChildren = rightNode->children.size();
        for( int c=0; c<numRightChildren; ++c )
            rightNode->children[c] = lastIndex - rightNode->children[c];
    }
    // Handle the middle node if it exists
    if( numLocalNodes % 2 != 0 )
    {
        const int midIndex = numLocalNodes/2;
        // Update the parent indices to the final values
        LocalSymmNode* middleNode = eTree.localNodes[midIndex];
        LocalSepOrLeaf* middleSepOrLeaf = sepTree.localSepsAndLeaves[midIndex];
        middleNode->parent = lastIndex - middleNode->parent;
        middleSepOrLeaf->parent = lastIndex - middleSepOrLeaf->parent;
        // Update the children's indices
        const int numChildren = middleNode->children.size();
        for( int c=0; c<numChildren; ++c )
            middleNode->children[c] = lastIndex - middleNode->children[c];
    }
}

inline void
BuildChildrenFromPerm
( const Graph& graph, const std::vector<int>& perm, 
  int leftChildSize, Graph& leftChild,
  int rightChildSize, Graph& rightChild )
{
#ifndef RELEASE
    PushCallStack("BuildChildrenFromPerm");
#endif
    const int numSources = graph.NumSources();
    const int sepSize = numSources - leftChildSize - rightChildSize;

    // Build the inverse permutation
    std::vector<int> inversePerm( numSources );
    for( int i=0; i<numSources; ++i )
        inversePerm[perm[i]] = i;

    // Get an upper bound on the number of edges in the child graphs
    int leftChildUpperBound=0, rightChildUpperBound=0;
    for( int s=0; s<leftChildSize; ++s )
        leftChildUpperBound += graph.NumConnections( inversePerm[s] );
    for( int s=0; s<rightChildSize; ++s )
        rightChildUpperBound += 
            graph.NumConnections( inversePerm[s+leftChildSize] );

    // Build the left child's graph
    leftChild.ResizeTo( leftChildSize );
    leftChild.StartAssembly();
    leftChild.Reserve( leftChildUpperBound );
    for( int s=0; s<leftChildSize; ++s )
    {
        const int source = s;
        const int inverseSource = inversePerm[s];
        const int offset = graph.EdgeOffset( inverseSource );
        const int numConnections = graph.NumConnections( inverseSource );
        for( int t=0; t<numConnections; ++t )
        {
            const int inverseTarget = graph.Target( offset+t );
            const int target = ( inverseTarget < numSources ? 
                                 perm[inverseTarget] :
                                 inverseTarget );
#ifndef RELEASE
            if( target >= leftChildSize && target < (numSources-sepSize) )
                throw std::logic_error
                ("Invalid bisection, left set touches right set");
#endif
            leftChild.Insert( source, target );
        }
    }
    leftChild.StopAssembly();

    // Build the right child's graph
    rightChild.ResizeTo( rightChildSize );
    rightChild.StartAssembly();
    rightChild.Reserve( rightChildUpperBound );
    for( int s=0; s<rightChildSize; ++s )
    {
        const int source = s+leftChildSize;
        const int inverseSource = inversePerm[source];
        const int offset = graph.EdgeOffset( inverseSource );
        const int numConnections = graph.NumConnections( inverseSource );
        for( int t=0; t<numConnections; ++t )
        {
            const int inverseTarget = graph.Target( offset+t );
            const int target = ( inverseTarget < numSources ?
                                 perm[inverseTarget] :
                                 inverseTarget );
#ifndef RELEASE
            if( target < leftChildSize )
                throw std::logic_error
                ("Invalid bisection, right set touches left set");
#endif
            // The targets that are in parent separators do not need to be
            rightChild.Insert( source-leftChildSize, target-leftChildSize );
        }
    }
    rightChild.StopAssembly();
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void 
BuildChildFromPerm
( const DistGraph& graph, const DistMap& perm,
  int leftChildSize, int rightChildSize,
  bool& onLeft, DistGraph& child )
{
#ifndef RELEASE
    PushCallStack("BuildChildFromPerm");
#endif
    const int numSources = graph.NumSources();
    const int numLocalSources = graph.NumLocalSources();
    const int sepSize = numSources - leftChildSize - rightChildSize;

    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    // Build the child graph from the partitioned parent
    const int smallTeamSize = commSize/2;
    const int largeTeamSize = commSize - smallTeamSize;
    const bool inSmallTeam = ( commRank < smallTeamSize );
    const bool smallOnLeft = ( leftChildSize <= rightChildSize );
    const int leftTeamSize = ( smallOnLeft ? smallTeamSize : largeTeamSize );
    const int rightTeamSize = ( smallOnLeft ? largeTeamSize : smallTeamSize );
    const int leftTeamOffset = ( smallOnLeft ? 0 : smallTeamSize );
    const int rightTeamOffset = ( smallOnLeft ? smallTeamSize : 0 );
    onLeft = ( inSmallTeam == smallOnLeft );

    const int leftTeamBlocksize = leftChildSize / leftTeamSize;
    const int rightTeamBlocksize = rightChildSize / rightTeamSize;

    // Count how many rows we must send to each process 
    std::vector<int> rowSendSizes( commSize, 0 );
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOffset + 
                RowToProcess( i, leftTeamBlocksize, leftTeamSize );
            ++rowSendSizes[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOffset +
                RowToProcess
                ( i-leftChildSize, rightTeamBlocksize, rightTeamSize );
            ++rowSendSizes[q];
        }
    }

    // Exchange the number of rows
    std::vector<int> rowRecvSizes( commSize );
    mpi::AllToAll( &rowSendSizes[0], 1, &rowRecvSizes[0], 1, comm );

    // Prepare for the AllToAll to exchange the row indices and 
    // the number of column indices per row
    int numSendRows=0;
    std::vector<int> rowSendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        rowSendOffsets[q] = numSendRows;
        numSendRows += rowSendSizes[q];
    }
    int numRecvRows=0;
    std::vector<int> rowRecvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        rowRecvOffsets[q] = numRecvRows;
        numRecvRows += rowRecvSizes[q];
    }

    // Pack the row indices and how many column entries there will be per row
    std::vector<int> rowSendLengths( numSendRows );
    std::vector<int> rowSendIndices( numSendRows );
    std::vector<int> offsets = rowSendOffsets;
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOffset + 
                RowToProcess( i, leftTeamBlocksize, leftTeamSize );
            rowSendIndices[offsets[q]] = i;
            rowSendLengths[offsets[q]] = graph.NumConnections( s );
            ++offsets[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOffset + 
                RowToProcess
                ( i-leftChildSize, rightTeamBlocksize, rightTeamSize );
            rowSendIndices[offsets[q]] = i;
            rowSendLengths[offsets[q]] = graph.NumConnections( s );
            ++offsets[q];
        }
    }

    // Perform the row lengths exchange
    std::vector<int> rowRecvLengths( numRecvRows );
    mpi::AllToAll
    ( &rowSendLengths[0], &rowSendSizes[0], &rowSendOffsets[0],
      &rowRecvLengths[0], &rowRecvSizes[0], &rowRecvOffsets[0], comm );

    // Perform the row indices exchange
    std::vector<int> rowRecvIndices( numRecvRows );
    mpi::AllToAll
    ( &rowSendIndices[0], &rowSendSizes[0], &rowSendOffsets[0],
      &rowRecvIndices[0], &rowRecvSizes[0], &rowRecvOffsets[0], comm );
    rowSendIndices.clear();

    // Set up for sending the column indices
    int numSendIndices=0;
    std::vector<int> indexSendSizes( commSize, 0 );
    std::vector<int> indexSendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowSendSizes[q];
        const int offset = rowSendOffsets[q];
        for( int s=0; s<numRows; ++s )
            indexSendSizes[q] += rowSendLengths[offset+s];

        indexSendOffsets[q] = numSendIndices;
        numSendIndices += indexSendSizes[q];
    }
    rowSendLengths.clear();
    int numRecvIndices=0;
    std::vector<int> indexRecvSizes( commSize, 0 );
    std::vector<int> indexRecvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowRecvSizes[q];
        const int offset = rowRecvOffsets[q];
        for( int s=0; s<numRows; ++s )
            indexRecvSizes[q] += rowRecvLengths[offset+s];

        indexRecvOffsets[q] = numRecvIndices;
        numRecvIndices += indexRecvSizes[q];
    }
    rowSendSizes.clear();
    rowSendOffsets.clear();

    // Pack the indices
    std::vector<int> sendIndices( numSendIndices );
    offsets = indexSendOffsets;
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOffset + 
                RowToProcess( i, leftTeamBlocksize, leftTeamSize );

            int& offset = offsets[q];
            const int numConnections = graph.NumConnections( s );
            const int localEdgeOffset = graph.LocalEdgeOffset( s );
            for( int j=0; j<numConnections; ++j )
                sendIndices[offset++] = graph.Target( localEdgeOffset+j );
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOffset + 
                RowToProcess
                ( i-leftChildSize, rightTeamBlocksize, rightTeamSize );
               
            int& offset = offsets[q];
            const int numConnections = graph.NumConnections( s );
            const int localEdgeOffset = graph.LocalEdgeOffset( s );
            for( int j=0; j<numConnections; ++j )
                sendIndices[offset++] = graph.Target( localEdgeOffset+j );
        }
    }

    // Send/recv the column indices
    std::vector<int> recvIndices( numRecvIndices );
    mpi::AllToAll
    ( &sendIndices[0], &indexSendSizes[0], &indexSendOffsets[0],
      &recvIndices[0], &indexRecvSizes[0], &indexRecvOffsets[0], comm );
    sendIndices.clear();
    indexSendSizes.clear();
    indexSendOffsets.clear();

    // Get the indices after reordering
    perm.Translate( recvIndices );

    // Put the connections into our new graph
    const int childTeamRank = 
        ( onLeft ? commRank-leftTeamOffset : commRank-rightTeamOffset );
    MPI_Comm childComm;
    mpi::CommSplit( comm, onLeft, childTeamRank, childComm );
    child.SetComm( childComm );
    if( onLeft )
        child.ResizeTo( leftChildSize );
    else
        child.ResizeTo( rightChildSize );

    child.StartAssembly();
    child.Reserve( recvIndices.size() );
    int offset=0;
    for( int s=0; s<numRecvRows; ++s )
    {
        const int source = rowRecvIndices[s];
        const int numConnections = rowRecvLengths[s];
#ifndef RELEASE
        const int childFirstLocalSource = child.FirstLocalSource();
        const int numChildLocalSources = child.NumLocalSources();
        const int adjustedSource = ( onLeft ? source : source-leftChildSize );
        if( adjustedSource < childFirstLocalSource || 
            adjustedSource >= childFirstLocalSource+numChildLocalSources )
            throw std::logic_error("source was out of bounds");
#endif
        for( int t=0; t<numConnections; ++t )
        {
            const int target = recvIndices[offset++];
            if( onLeft )
            {
#ifndef RELEASE
                if( target >= leftChildSize && target < (numSources-sepSize) )
                {
                    std::ostringstream msg;
                    msg << "Invalid dist bisection, left set touches right:\n"
                        << "  " << source << " touches " << target << " and "
                        << "leftChildSize=" << leftChildSize;
                    throw std::logic_error( msg.str().c_str() );
                }
#endif
                child.Insert( source, target );
            }
            else
            {
#ifndef RELEASE
                if( target < leftChildSize )
                    throw std::logic_error
                    ("Invalid dist bisection, right set touches left set");
#endif
                child.Insert( source-leftChildSize, target-leftChildSize );
            }
        }
    }
    child.StopAssembly();
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
BuildMap
( const DistGraph& graph, 
  const DistSeparatorTree& sepTree, 
        DistMap& map )
{
#ifndef RELEASE
    PushCallStack("BuildMap");
#endif
    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int numSources = graph.NumSources();

    map.SetComm( comm );
    map.ResizeTo( numSources );
    const int blocksize = map.Blocksize();

    const int numLocal = sepTree.localSepsAndLeaves.size();
    // NOTE: The dist separator tree does not double-count the first 
    //       single-process node, but DistSymmInfo does. Thus their number of
    //       distributed nodes is different by one.
    const int numDist = sepTree.distSeps.size();

    // Traverse local sepTree to count how many indices we should send the
    // final index for
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const int numIndices = sepOrLeaf.indices.size();
        for( int t=0; t<numIndices; ++t )
        {
            const int i = sepOrLeaf.indices[t];
#ifndef RELEASE
            if( i < 0 || i >= graph.NumSources() )
            {
                std::ostringstream msg;
                msg << "local separator index, " << i << ", was not in [0,"
                    << graph.NumSources() << ")";
                throw std::logic_error( msg.str().c_str() );
            }
#endif
            const int q = RowToProcess( i, blocksize, commSize );
            ++sendSizes[q];
        }
    }
    for( int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const int numIndices = sep.indices.size();
        const int teamSize = mpi::CommSize( sep.comm );
        const int teamRank = mpi::CommRank( sep.comm );
        const int numLocalIndices = 
            LocalLength( numIndices, teamRank, teamSize );
        for( int tLocal=0; tLocal<numLocalIndices; ++tLocal )
        {
            const int t = teamRank + tLocal*teamSize;
            const int i = sep.indices[t];
#ifndef RELEASE
            if( i < 0 || i >= graph.NumSources() )
            {
                std::ostringstream msg;
                msg << "dist separator index, " << i << ", was not in [0,"
                    << graph.NumSources() << ")";
                throw std::logic_error( msg.str().c_str() );
            }
#endif
            const int q = RowToProcess( i, blocksize, commSize );
            ++sendSizes[q];
        }
    }

    // Use a single-entry AllToAll to coordinate how many indices will be 
    // exchanges
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm );

    // Pack the reordered indices
    int numSends = 0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSends;
        numSends += sendSizes[q];
    }
    std::vector<int> sendIndices( numSends );
    std::vector<int> sendOrigIndices( numSends );
    std::vector<int> offsets = sendOffsets;
    for( int s=0; s<numLocal; ++s )
    {
        const LocalSepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const int numIndices = sepOrLeaf.indices.size();
        for( int t=0; t<numIndices; ++t )
        {
            const int i = sepOrLeaf.indices[t];
            const int iMapped = sepOrLeaf.offset + t;
            const int q = RowToProcess( i, blocksize, commSize );
            sendOrigIndices[offsets[q]] = i;
            sendIndices[offsets[q]] = iMapped;
            ++offsets[q];
        }
    }
    for( int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const int numIndices = sep.indices.size();
        const int teamSize = mpi::CommSize( sep.comm );
        const int teamRank = mpi::CommRank( sep.comm );
        const int numLocalIndices = 
            LocalLength( numIndices, teamRank, teamSize );
        for( int tLocal=0; tLocal<numLocalIndices; ++tLocal )
        {
            const int t = teamRank + tLocal*teamSize;
            const int i = sep.indices[t];
            const int iMapped = sep.offset + t;
            const int q = RowToProcess( i, blocksize, commSize );
            sendOrigIndices[offsets[q]] = i;
            sendIndices[offsets[q]] = iMapped;
            ++offsets[q];
        }
    }

    // Perform an AllToAll to exchange the reordered indices
    int numRecvs = 0;
    std::vector<int> recvOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvOffsets[q] = numRecvs;
        numRecvs += recvSizes[q];
    }
#ifndef RELEASE
    const int numLocalSources = map.NumLocalSources();
    if( numRecvs != numLocalSources )
        throw std::logic_error("incorrect number of recv indices");
#endif
    std::vector<int> recvIndices( numRecvs );
    mpi::AllToAll
    ( &sendIndices[0], &sendSizes[0], &sendOffsets[0],
      &recvIndices[0], &recvSizes[0], &recvOffsets[0], comm );

    // Perform an AllToAll to exchange the original indices
    std::vector<int> recvOrigIndices( numRecvs );
    mpi::AllToAll
    ( &sendOrigIndices[0], &sendSizes[0], &sendOffsets[0],
      &recvOrigIndices[0], &recvSizes[0], &recvOffsets[0], comm );

    // Unpack the indices
    const int firstLocalSource = graph.FirstLocalSource();
    for( int s=0; s<numRecvs; ++s )
    {
        const int i = recvOrigIndices[s];
        const int iLocal = i - firstLocalSource;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= numLocalSources )
        {
            std::ostringstream msg;
            msg << "Local index was out of bounds: " << iLocal 
                << " is not in [0," << numLocalSources << ")" << std::endl;
            throw std::logic_error( msg.str().c_str() );
        }
#endif
        const int iMapped = recvIndices[s];
#ifndef RELEASE
        if( iMapped < 0 || iMapped >= graph.NumSources() )
        {
            std::ostringstream msg;
            msg << "mapped index, " << iMapped << ", was not in [0,"
                << graph.NumSources() << ")";
            throw std::logic_error( msg.str().c_str() );
        }
#endif
        map.SetLocal( iLocal, iMapped );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
