Nested dissection
=================
*Nested dissection* refers to the process of recursively partitioning a graph
into two disjoint subdomains which only interact through a region known as a 
*separator* (for :math:`d`-dimensional domains, separators are typically 
:math:`d-1`-dimensional); each separator from the nested dissection process 
yields a node in a *separator tree*, which can be used to guide a multifrontal 
algorithm. The following routine uses a parallel graph partitioner (ParMETIS)
as a means of producing such a separator tree from an arbitrary graph. 

.. cpp:function:: void NestedDissection( const DistGraph& graph, std::vector<int>& localMap, DistSeparatorTree& sepTree, DistSymmInfo& info, int cutoff=128, int numDistSeps=10, int numSeqSeps=5, bool storeFactRecvIndices=true )

   In addition to generating the separator tree, which simply stores the indices
   of each separator, each process computes its portion of the mapping from the 
   original to the reordered indices (implied by nested dissection) within the 
   ``localMap`` vector, and the symbolic factorization, as well as a variety of 
   other useful information, is returned in an instance of the ``DistSymmInfo``
   structure, ``info``. The ``cutoff`` parameter determines the maximum 
   acceptable leaf node size for nested dissection, the ``numDistSeps`` and 
   ``numSeqSeps`` variables respectively determine how many distributed and 
   sequential separators should be tried for each bisection, and 
   ``storeFactRecvIndices`` determines whether or not to store information
   needed for the redistributions which occur in the subsequent numerical 
   factorization

   See `tests/NestedDissection <https://github.com/poulson/Clique/blob/master/tests/NestedDissection.cpp>`__ for an example of its usage directly on a 
   distributed graph, and `tests/DistSparseMatrix.cpp <https://github.com/poulson/Clique/blob/master/tests/DistSparseMatrix.cpp>`__ for its application to the
   underlying graph of a sparse matrix.
