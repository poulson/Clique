Factorization-related data structures
=====================================
**TODO**

DistSymmFrontTree
-----------------

.. cpp:type:: enum SolveMode

   Can be set to either
   
   * ``NORMAL_1D``

   * ``NORMAL_2D``

   * ``FAST_1D_LDL``

   * ``FAST_2D_LDL``

.. cpp:type:: struct LocalSymmFront<F>

   **TODO: Briefly describe each member variable**

   .. cpp:member:: Matrix<F> frontL

   .. cpp:member:: mutable Matrix<F> work

.. cpp:type:: struct DistSymmFront<F>

   **TODO: Briefly describe each member variable**

   .. cpp:member:: mutable DistMatrix<F,VC,STAR> front1dL

   .. cpp:member:: mutable DistMatrix<F,VC,STAR> work1d

   .. cpp:member:: mutable DistMatrix<F> front2dL

   .. cpp:member:: mutable DistMatrix<F> work2d

   .. cpp:member:: DistMatrix<F,VC,STAR> diag

.. cpp:type:: DistSymmFrontTree<F>

   **TODO: Briefly describe each member variable**

   .. cpp:member:: SolveMode mode

   .. cpp:member:: std::vector<LocalSymmFront<F> > localFronts

   .. cpp:member:: std::vector<DistSymmFront<F> > distFronts

   .. cpp:function:: DistSymmFrontTree( Orientation orientation, const DistSparseMatrix<F>& A, const std::vector<int>& localMap, const DistSeparatorTree& sepTree, const DistSymmInfo& info )
