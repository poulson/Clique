Analysis
========
Once the elimination tree has been decided, it is necessary to determine
the fill-in that will occur during factorization; that is to say, we must
determine which degrees of freedom each (super)node of the elimination tree is
coupled to after factorization. The set of vertices that a (super)node is
coupled to is commonly referred to as its *structure*, and
*symbolic factorization* consists of converting the original structure of each
node in the elimination tree of a sparse matrix into the resulting structure
after factorization.

The following routine is called ``SymmetricAnalysis`` instead of
``SymmetricSymbolicFactorization`` because, in addition to performing the
symbolic factorization, it also computes and stores other data which is
useful for numerical factorization and solves.

.. cpp:function:: void SymmetricAnalysis( const DistSymmElimTree& eTree, DistSymmInfo& info, bool storeFactRecvIndices=true )
    
.. note:: 
   Most users will not need to directly call this routine, as it is 
   part of :cpp:func:`NestedDissection`. It will only need to be manually 
   called in cases where someone has manually computed the elimination tree of 
   their sparse matrix.

Data structures
---------------
**The data structures used for symbolic analysis are meant to serve as 
black-boxes for the vast majority of users.** Nevertheless, they are 
briefly documented here for the sake of posterity.

DistSymmElimTree
^^^^^^^^^^^^^^^^

.. cpp:type:: struct SymmNode

   .. cpp:member:: int size

      The size of this node in the elimination tree.

   .. cpp:member:: int offset

      The first reordered index of the vertices in this node.

   .. cpp:member:: int parent

      The index of the local parent (:math:`-1` if it does not exist).

   .. cpp:member:: std::vector<int> children

      The indices of the child nodes.

   .. cpp:member:: std::vector<int> lowerStruct

      The sorted reordered indices of the connections to ancestor nodes.

.. cpp:type:: struct DistSymmNode

   .. cpp:member:: bool onLeft

      Whether or not this node is the left child of its parent 
      (assuming it has a parent).

   .. cpp:member:: mpi::Comm comm

      The communicator for the distributed node.

   .. cpp:member:: int size

      The number of vertices in this node.

   .. cpp:member:: int offset

      The first reordered index for the vertices in this node.

   .. cpp:member:: std::vector<int> lowerStruct

      The sorted reordered indices of the connections to ancestor nodes.

.. cpp:type:: struct DistSymmElimTree

   .. cpp:member:: std::vector<SymmNode*> localNodes

      An array of *pointers* to the local nodes.

   .. cpp:member:: std::vector<DistSymmNode> distNodes

      An array of distributed nodes, including the single-process 
      separator or leaf.

