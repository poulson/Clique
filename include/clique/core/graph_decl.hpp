/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

class Graph
{
public:
    // Constructors and destructors
    Graph();
    Graph( int numVertices );
    Graph( int numSources, int numTargets );
    Graph( const Graph& graph );
    // NOTE: This requires the DistGraph to be over a single process
    Graph( const DistGraph& graph );
    ~Graph();

    // High-level information
    int NumSources() const;
    int NumTargets() const;

    // Assembly-related routines
    void StartAssembly();
    void StopAssembly();
    void Reserve( int numEdges );
    void Insert( int source, int target );
    int Capacity() const;

    // Data
    int NumEdges() const;
    int Source( int edge ) const;
    int Target( int edge ) const;
    int EdgeOffset( int source ) const;
    int NumConnections( int source ) const;

    // For resizing the graph
    void Empty();
    void ResizeTo( int numVertices );
    void ResizeTo( int numSources, int numTargets );

    // For copying one graph into another
    const Graph& operator=( const Graph& graph );
    // NOTE: This requires the DistGraph to be over a single process
    const Graph& operator=( const DistGraph& graph );

    void Print( std::string msg ) const;

private:
    int numSources_, numTargets_;
    std::vector<int> sources_, targets_;

    // Helpers for local indexing
    bool assembling_, sorted_;
    std::vector<int> edgeOffsets_;
    void ComputeEdgeOffsets();

    static bool ComparePairs
    ( const std::pair<int,int>& a, const std::pair<int,int>& b );

    void EnsureNotAssembling() const;
    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    friend class DistGraph;
    template<typename F> friend class SparseMatrix;
};

} // namespace cliq
