Nested dissection
=================
*Nested dissection* refers to the process of recursively partitioning a graph
into two disjoint subdomains which only interact through a region known as a
*separator* (for :math:`d`-dimensional domains, separators are typically
:math:`d-1`-dimensional); each separator from the nested dissection process
yields a node in a *separator tree*, which can be used to guide a multifrontal
algorithm. The following routine uses a parallel graph partitioner (ParMETIS)
as a means of producing such a separator tree from an arbitrary graph.

.. cpp:function:: void NestedDissection( const DistGraph& graph, DistMap& map, DistSeparatorTree& sepTree, DistSymmInfo& info, bool sequential=true, int cutoff=128, int numDistSeps=1, int numSeqSeps=1, bool storeFactRecvInd=true )

   .. note:: 

      Since partition refinement via KL-FM is not very scalable from the runtime
      perspective, by default, graph bisection happens on a single processor 
      for each node of the frontal tree.

   In addition to generating the separator tree, which simply stores the indices
   of each separator, each process computes its portion of the mapping from the
   original to the reordered indices (implied by nested dissection) within the
   `map` vector, and the symbolic factorization, as well as a variety of
   other useful information, is returned in an instance of the 
   :cpp:type:`DistSymmInfo` structure, `info`. 
   The `cutoff` parameter determines the maximum
   acceptable leaf node size for nested dissection, the `numDistSeps` and
   `numSeqSeps` variables respectively determine how many distributed and
   sequential separators should be tried for each bisection, and
   `storeFactRecvInd` determines whether or not to store information
   needed for the redistributions which occur in the subsequent numerical
   factorization

   See `tests/NestedDissection <https://github.com/poulson/Clique/blob/master/tests/NestedDissection.cpp>`__ for an example of its usage directly on a
   distributed graph, and `tests/Solve.cpp <https://github.com/poulson/Clique/blob/master/tests/Solve.cpp>`__ for its application to the
   underlying graph of a sparse matrix.

.. cpp:function:: void NaturalNestedDissection( int nx, int ny, int nz, const DistGraph& graph, DistMap& map, DistSeparatorTree& sepTree, DistSymmInfo& info, int cutoff=128, bool storeFactRecvInd=true )

   Similar to :cpp:func:`NestedDissection`, but this version is specialized for 
   regular 3D grids where vertices are only connected to their nearest 
   neighbors. In this case, the graph can analytically be recursively bisected,
   and so the difficulties in parallelizing the KL-FM refinement can be avoided.

   See `tests/NaturalSolve <https://github.com/poulson/Clique/blob/master/tests/NaturalSolve.cpp>`__ for an example.

Data structures
---------------
**The data structures used for nested dissection are meant to serve as 
black-boxes for the vast majority of users.** Nevertheless, they are 
briefly documented here for the sake of posterity.

DistSeparatorTree
^^^^^^^^^^^^^^^^^

.. cpp:type:: struct SepOrLeaf

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

   .. cpp:member:: std::vector<SepOrLeaf*> localSepsAndLeaves

      An array of *pointers* to local separators and leaves.

   .. cpp:member:: std::vector<DistSeparator> distSeps

      An array of distributed separators.

      .. note::

         This array does *not* include the single process separator/leaf.

DistSymmInfo
^^^^^^^^^^^^

.. cpp:type:: struct SymmNodeInfo

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

   .. cpp:member:: bool onLeft

      Whether or not this node is a left child (assuming it has a parent).

   .. cpp:member:: int myOffset

      The sum of the node sizes for all previously ordered nodes.

   .. cpp:member:: std::vector<int> lowerStruct

      The sorted reordered indices of this node's connections to its ancestors
      **after factorization**.

   .. cpp:member:: std::vector<int> origLowerRelInd

      Maps from the original lower structure to their placement in the 
      structure after factorization.

   .. cpp:member:: std::vector<int> leftRelInd
   .. cpp:member:: std::vector<int> rightRelInd

      The relative indices of the left/right child's lower structure into this 
      structure.

.. cpp:type:: struct FactorMetadata

   .. cpp:member:: std::vector<int> numChildSendInd

   .. cpp:member:: std::deque<int> leftColInd
   .. cpp:member:: std::deque<int> leftRowInd
   .. cpp:member:: std::deque<int> rightColInd
   .. cpp:member:: std::deque<int> rightRowInd

   .. cpp:member:: mutable std::vector<std::deque<int> > childRecvInd

   .. cpp:function:: void EmptyChildRecvIndices() const

      Clears ``childRecvInd``

   .. cpp:function:: void Empty()

      Clears all members of structure

.. cpp:type:: struct SolveMetadata1d

   .. cpp:member:: int localSize

   .. cpp:member:: std::deque<int> leftInd
   .. cpp:member:: std::deque<int> rightInd

   .. cpp:member:: std::vector<int> numChildSendInd

   .. cpp:member:: std::vector<std::deque<int> > childRecvInd

   .. cpp:function:: void Empty()

      Clears all members of structure

.. cpp:type:: struct SolveMetadata2d

   .. cpp:member:: int localHeight
   .. cpp:member:: int localWidth

   .. cpp:member:: std::deque<int> leftInd
   .. cpp:member:: std::deque<int> rightInd

   .. cpp:member:: std::vector<int> numChildSendInd

   .. cpp:member:: std::vector<std::deque<int> > childRecvInd

   .. cpp:function:: void Empty()

      Clears all members of structure

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

   .. cpp:member:: int leftSize
   .. cpp:member:: int rightSize

      The number of vertices in the left/right child (assuming it exists).

   .. cpp:member:: std::vector<int> lowerStruct

      The sorted reordered indices of this node's connections to its 
      ancestors **after factorization**.

   .. cpp:member:: std::vector<int> origLowerRelInd

      Maps from the original lower structure to their placement in the 
      structure after factorization.

   .. cpp:member:: std::vector<int> leftRelInd
   .. cpp:member:: std::vector<int> rightRelInd

      The relative indices of the left/right child's lower structure into this 
      structure.

   .. cpp:member:: FactorMetadata factorMeta
   .. cpp:member:: SolveMetadata1d solveMeta1d
   .. cpp:member:: SolveMetadata2d solveMeta2d

.. cpp:type:: struct DistSymmInfo

   .. cpp:member:: std::vector<SymmNodeInfo> localNodes

   .. cpp:member:: std::vector<DistSymmNodeInfo> distNodes

