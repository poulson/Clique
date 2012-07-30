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
   `localMap` vector, and the symbolic factorization, as well as a variety of
   other useful information, is returned in an instance of the 
   :cpp:type:`DistSymmInfo` structure, `info`. 
   The `cutoff` parameter determines the maximum
   acceptable leaf node size for nested dissection, the `numDistSeps` and
   `numSeqSeps` variables respectively determine how many distributed and
   sequential separators should be tried for each bisection, and
   `storeFactRecvIndices` determines whether or not to store information
   needed for the redistributions which occur in the subsequent numerical
   factorization

   See `tests/NestedDissection <https://github.com/poulson/Clique/blob/master/tests/NestedDissection.cpp>`__ for an example of its usage directly on a
   distributed graph, and `tests/VectorSolve.cpp <https://github.com/poulson/Clique/blob/master/tests/VectorSolve.cpp>`__ for its application to the
   underlying graph of a sparse matrix.

Data structures
---------------
**The data structures used for nested dissection are meant to serve as 
black-boxes for the vast majority of users.** Nevertheless, they are 
briefly documented here for the sake of posterity.

DistSeparatorTree
^^^^^^^^^^^^^^^^^

.. cpp:type:: struct LocalSepOrLeaf

   .. cpp:member:: int parent 

      The index of the local parent separator (:math:`-1` if it does not exist).

   .. cpp:member:: int offset

      The first reordered index of this node.

   .. cpp:member:: std::vector<int> indices

      The original indices of this node.

.. cpp:type:: struct DistSeparator

   .. cpp:member:: mpi::Comm comm

      The communicator for this distributed separator.

   .. cpp:member:: int offset

      The first reordered index of this separator.

   .. cpp:member:: std::vector<int> indices

      The original indices of this separator.

.. cpp:type:: struct DistSeparatorTree

   .. cpp:member:: std::vector<LocalSepOrLeaf*> localSepsAndLeaves

      An array of *pointers* to local separators and leaves.

   .. cpp:member:: std::vector<DistSeparator> distSeps

      An array of distributed separators.

      .. note::

         This array does *not* include the single process separator/leaf.

DistSymmInfo
^^^^^^^^^^^^

.. cpp:type:: struct LocalSymmNodeInfo

   .. rubric:: Known before analysis

   .. cpp:member:: int size

      The number of vertices in this node.

   .. cpp:member:: int offset

      The first reordered index of the vertices in this node.

   .. cpp:member:: std::vector<int> children

      The indices of the child nodes.

   .. cpp:member:: std::vector<int> origLowerStruct

      The original sorted reordered indices of this node's connections to its
      ancestors.

   .. rubric:: Computed during analysis

   .. cpp:member:: bool isLeftChild

      Whether or not this node is a left child (assuming it has a parent).

   .. cpp:member:: int myOffset

      The sum of the node sizes for all previously ordered nodes.

   .. cpp:member:: std::vector<int> lowerStruct

      The sorted reordered indices of this node's connections to its ancestors
      **after factorization**.

   .. cpp:member:: std::vector<int> origLowerRelIndices

      Maps from the original lower structure to their placement in the 
      structure after factorization.

   .. cpp:member:: std::vector<int> leftChildRelIndices

      The relative indices of the left child's lower structure into this 
      structure.

   .. cpp:member:: std::vector<int> rightChildRelIndices

      The relative indices of the right child's lower structure into this 
      structure.

.. cpp:type:: struct DistSymmNodeInfo

   .. rubric:: Known before analysis

   .. cpp:member:: int size

      The number of vertices in this node.

   .. cpp:member:: int offset

      The first reordered index of the vertices in this node.

   .. cpp:member:: std::vector<int> origLowerStruct

      The original sorted reordered indices of this node's connections to its
      ancestors.

   .. cpp:member:: bool onLeft

      Whether or not this node is a left child (assuming it has a parent).

   .. cpp:member:: mpi::Comm comm

      The communicator for this leaf or separator.

   .. rubric:: Computed during analysis

   .. cpp:member:: Grid* grid

      The process grid which will be used to distribute the frontal matrix for
      this node.

   .. cpp:member:: int myOffset

      The sum of the node sizes for all previously ordered nodes.

   .. cpp:member:: int leftChildSize

      The number of vertices in the left child (assuming it exists).

   .. cpp:member:: int rightChildSize

      The number of vertices in the right child (assuming it exists).

   .. cpp:member:: std::vector<int> lowerStruct

      The sorted reordered indices of this node's connections to its 
      ancestors **after factorization**.

   .. cpp:member:: std::vector<int> origLowerRelIndices

      Maps from the original lower structure to their placement in the 
      structure after factorization.

   .. cpp:member:: std::vector<int> leftChildRelIndices

      The relative indices of the left child's lower structure into this 
      structure.

   .. cpp:member:: std::vector<int> rightChildRelIndices

      The relative indices of the right child's lower structure into this 
      structure.

   .. cpp:member:: std::vector<int> numChildFactSendIndices

      **Left off here**

   .. cpp:member:: std::vector<int> leftChildFactColIndices

   .. cpp:member:: std::vector<int> leftChildFactRowIndices

   .. cpp:member:: std::vector<int> rightChildFactColIndices

   .. cpp:member:: std::vector<int> rightChildFactRowIndices

   .. cpp:member:: mutable std::vector<std::deque<int> > childFactRecvIndices

   .. cpp:member:: std::deque<int> leftChildSolveIndices

   .. cpp:member:: std::deque<int> rightChildSolveIndices

   .. cpp:member:: int localSize1d

   .. cpp:member:: int localOffset1d

   .. cpp:member:: std::vector<int> numChildSolveSendIndices

   .. cpp:member:: std::vector<std::deque<int> > childSolveRecvIndices

.. cpp:type:: struct DistSymmInfo

   .. cpp:member:: std::vector<LocalSymmNodeInfo> localNodes

   .. cpp:member:: std::vector<DistSymmNodeInfo> distNodes

