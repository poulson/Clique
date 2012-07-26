Analysis-related data structures
================================
**The data structures used for symbolic analysis are meant to serve as 
black-boxes for the vast majority of users.** Nevertheless, they are 
briefly documented here for the sake of posterity.

DistSeparatorTree
-----------------

.. cpp:type:: struct LocalSepOrLeaf

   **TODO: Briefly describe each member variable**

   .. cpp:member:: int parent

   .. cpp:member:: int offset

   .. cpp:member:: std::vector<int> indices

.. cpp:type:: struct DistSeparator

   .. cpp:member:: mpi::Comm comm

   .. cpp:member:: int offset

   .. cpp:member:: std::vector<int> indices

.. cpp:type:: struct DistSeparatorTree

   .. cpp:member:: std::vector<LocalSepOrLeaf*> localSepsAndLeaves

   .. cpp:member:: std::vector<DistSeparator> distSeps

DistSymmElimTree
----------------

.. cpp:type:: struct LocalSymmNode

   **TODO: Briefly describe each member variable**

   .. cpp:member:: int size

   .. cpp:member:: int offset

   .. cpp:member:: int parent

   .. cpp:member:: std::vector<int> children

   .. cpp:member:: std::vector<int> lowerStruct

.. cpp:type:: struct DistSymmNode

   .. cpp:member:: bool onLeft

   .. cpp:member:: mpi::Comm comm

   .. cpp:member:: int size

   .. cpp:member:: int offset

   .. cpp:member:: std::vector<int> lowerStruct

.. cpp:type:: struct DistSymmElimTree

   .. cpp:member:: std::vector<LocalSymmNode*> localNodes

   .. cpp:member:: std::vector<DistSymmNode> distNodes

DistSymmInfo
------------

.. cpp:type:: struct LocalSymmNodeInfo

   **TODO: Briefly describe each member variable**

   .. rubric:: Known before analysis

   .. cpp:member:: int size

   .. cpp:member:: int offset

   .. cpp:member:: std::vector<int> children

   .. cpp:member:: std::vector<int> origLowerStruct

   .. rubric:: Computed during analysis

   .. cpp:member:: bool isLefChild

   .. cpp:member:: int myOffset

   .. cpp:member:: std::vector<int> lowerStruct

   .. cpp:member:: std::vector<int> origLowerRelIndices

   .. cpp:member:: std::vector<int> leftChildRelIndices

   .. cpp:member:: std::vector<int> rightChildRelIndices

.. cpp:type:: struct DistSymmNodeInfo

   .. rubric:: Known before analysis

   .. cpp:member:: int size

   .. cpp:member:: int offset

   .. cpp:member:: std::vector<int> origLowerStruct

   .. cpp:member:: bool onLeft

   .. cpp:member:: mpi::Comm comm

   .. rubric:: Computed during analysis

   .. cpp:member:: Grid* grid

   .. cpp:member:: int myOffset

   .. cpp:member:: int leftChildSize

   .. cpp:member:: int rightChildSize

   .. cpp:member:: std::vector<int> lowerStruct

   .. cpp:member:: std::vector<int> origLowerRelIndices

   .. cpp:member:: std::vector<int> leftChildRelIndices

   .. cpp:member:: std::vector<int> rightChildRelIndices

   .. cpp:member:: std::vector<int> numChildFactSendIndices

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


