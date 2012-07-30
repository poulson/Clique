Factorization
=============

LDL factorization
-----------------

.. cpp:function:: void LDL( Orientation orientation, DistSymmInfo& info, DistSymmFrontTree<F>& L )

   Performs either an :math:`LDL^T` or :math:`LDL^H` factorization, depending 
   on whether `orientation` is set to ``TRANSPOSE`` or ``ADJOINT``. See 
   `tests/VectorSolve <https://github.com/poulson/Clique/blob/master/tests/VectorSolve.cpp>`__ for an example usage.

   .. note::

      This routine does not pivot, so it should be used with caution on matrices
      which are not Hermitian positive-definite.

.. cpp:function:: void BlockLDL( Orientation orientation, DistSymmInfo& info, DistSymmFrontTree<F>& L )

   Performs a block-diagonal :math:`LDL^T` or :math:`LDL^H` factorization, 
   depending on whether `orientation` is set to ``TRANSPOSE`` or ``ADJOINT``.

   .. note::

      This routine does not pivot, so it should be used with caution on matrices
      which are not Hermitian positive-definite.

Data structures
---------------

.. note::

   The following data structure is meant to be treated as a black-box by the 
   vast majority of users, but its member variables are documented here for 
   users who are interested.

DistSymmFrontTree
^^^^^^^^^^^^^^^^^

This data structure represents a distributed symmetric :math:`LDL^T` or 
:math:`LDL^H` factorization.

.. cpp:type:: struct LocalSymmFront<T>

   .. cpp:member:: Matrix<T> frontL

      The left half of a local symmetric frontal matrix. If `frontL` is 
      :math:`m \times n`, then the associated node has :math:`n` vertices 
      and the lower structure is of size :math:`m-n`.

   .. cpp:member:: mutable Matrix<T> work

      Used for temporary results during operations with this frontal matrix.

.. cpp:type:: struct DistSymmFront<T>

   .. cpp:member:: mutable DistMatrix<T,VC,STAR> front1dL

      The left half of a frontal matrix in a 1D distribution. If `front1dL` is 
      :math:`m \times n`, then the associated node has :math:`n` vertices 
      and the lower structure is of size :math:`m-n`.

   .. cpp:member:: mutable DistMatrix<T,VC,STAR> work1d

      Used for temporary results which are in a 1D distribution.

   .. cpp:member:: mutable DistMatrix<T> front2dL

      The left half of a frontal matrix in a 2D distribution. 

   .. cpp:member:: mutable DistMatrix<T> work2d

      Used for temporary results which are in a 2D distribution.

   .. cpp:member:: DistMatrix<T,VC,STAR> diag

      Used for storing the diagonal of the frontal matrix.

.. cpp:type:: enum SolveMode

   Can be set to either
   
   * ``NORMAL_1D``: frontal matrices are distributed in a 1D manner

   * ``NORMAL_2D``: frontal matrices are distributed in a 2D manner (default)

   * ``FAST_1D_LDL``: frontal matrices are inverted and distributed in a 1D 
     manner

   * ``FAST_2D_LDL``: frontal matrices are inverted and distributed in a 2D 
     manner

.. cpp:type:: struct DistSymmFrontTree<T>

   .. cpp:member:: SolveMode mode

      Determines which format the distributed frontal matrices are in.

   .. cpp:member:: std::vector<LocalSymmFront<T> > localFronts

      The vector of local frontal matrices.

   .. cpp:member:: std::vector<DistSymmFront<T> > distFronts

      The vector of distributed frontal matrices.

   .. cpp:function:: DistSymmFrontTree( Orientation orientation, const DistSparseMatrix<T>& A, const DistMap& map, const DistSeparatorTree& sepTree, const DistSymmInfo& info )

      A constructor which converts a distributed sparse matrix into a symmetric
      frontal tree which is ready for factorization (e.g., with 
      :cpp:func:`LDL` or :cpp:func:`BlockLDL`).

.. cpp:type:: struct DistSymmFrontTree<F>

   Same as above, but this implies that the underlying datatype `F` is a field.
