Analysis
========
**TODO: Describe the data structures and routines**

Nested dissection
=================
*Nested dissection* refers to the process of recursively partitioning a graph
into two disjoint subdomains which only interact through a region known as a 
*separator* (for :math:`d`-dimensional domains, separators are typically 
:math:`d-1`-dimensional). **TODO: Finish describing nested dissection**

.. cpp:function:: void NestedDissection( const DistGraph& graph, std::vector<int>& localMap, DistSeparatorTree& sepTree, DistSymmInfo& info, int cutoff=128, int numDistSeps=10, int numSeqSeps=5, bool storeFactRecvIndices=true )
