Factorization
=============

LDL factorization
-----------------

.. cpp:function:: void LDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType )

   Performs either an :math:`LDL^T` or :math:`LDL^H` factorization, depending 
   on whether `L` is marked as Hermitian. See 
   `tests/Solve <https://github.com/poulson/Clique/blob/master/tests/Solve.cpp>`__ for an example usage.

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

.. cpp:type:: struct SymmFront<T>

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

.. cpp:type:: enum SymmFrontType

   Can be set to either
   
   * ``SYMM_1D``: Symmetric/Hermitian fronts distributed in a 1D manner

   * ``SYMM_2D``: Symmetric/Hermitian fronts distributed in a 2D manner 

   * ``LDL_1D``: LDL factorization distributed in a 1D manner

   * ``LDL_2D``: LDL factorization distributed in a 2D manner 

   * ``LDL_SELINV_1D``: LDL factorization with inverted diagonal blocks 
     distributed in a 1D manner

   * ``LDL_SELINV_2D``: LDL factorization with inverted diagonal blocks 
     distributed in a 2D manner

   * ``BLOCK_LDL_2D``: Block LDL factorization with fronts distributed in a 
     2D manner

.. cpp:type:: struct DistSymmFrontTree<T>

   .. cpp:member:: bool isHermitian

      If true, the matrix is assumed to be Hermitian; otherwise, it is 
      treated as symmetric.

   .. cpp:member:: SymmFrontType frontType

      Specifies the form of the frontal matrices.

   .. cpp:member:: std::vector<SymmFront<T> > localFronts

      The vector of local frontal matrices.

   .. cpp:member:: std::vector<DistSymmFront<T> > distFronts

      The vector of distributed frontal matrices.

   .. cpp:function:: DistSymmFrontTree( Orientation orientation, const DistSparseMatrix<T>& A, const DistMap& map, const DistSeparatorTree& sepTree, const DistSymmInfo& info )

      A constructor which converts a distributed sparse matrix into a symmetric
      frontal tree which is ready for factorization (e.g., with 
      :cpp:func:`LDL` or :cpp:func:`BlockLDL`).

   .. cpp:function:: void Initialize( Orientation orientation, const DistSparseMatrix<T>& A, const DistMap& map, const DistSeparatorTree& sepTree, const DistSymmInfo& info )

      The same as the :cpp:func:`DistSymmFrontTree\<T>::DistSymmFrontTree`
      constructor, but callable after construction.

.. cpp:type:: struct DistSymmFrontTree<F>

   Same as above, but this implies that the underlying datatype `F` is a field.
