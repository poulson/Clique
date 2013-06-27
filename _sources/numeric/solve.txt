Solves
======

Black-box solvers
-----------------

.. cpp:function:: void SymmetricSolve( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )

   Overwrites :math:`x` with :math:`A^{-1} x`, where :math:`A` is assumed to be 
   symmetric. The main optional argument is 
   whether or not graph partitioning should be performed sequentially or in 
   parallel. For modest-to-large numbers of processes, one should only use 
   parallel graph partitioning if the entire graph for the sparse matrix cannot
   be stored on a single process. The remaining optional arguments determine 
   how many distributed and sequential separators should be tested for each 
   separator (`numDistSeps` and `numSeqSeps`) and what the maximum allowed 
   subdomain size is (`cutoff`). See
   `tests/SimpleSolve <https://github.com/poulson/Clique/blob/master/tests/SimpleSolve.cpp>`__ for an example usage.

.. cpp:function:: void HermitianSolve( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )

   Overwrites :math:`x` with :math:`A^{-1} x`, where :math:`A` is assumed to be
   Hermitian. The optional arguments are identical to those of 
   :cpp:func:`SymmetricSolve`.

Solving after factorization
---------------------------

.. cpp:function:: void Solve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )

   After having performed an :math:`LDL^T` or :math:`LDL^H` factorization 
   (via :cpp:func:`LDL`), this routine can be used to solve a set of 
   right-hand sides. See
   `tests/Solve <https://github.com/poulson/Clique/blob/master/tests/Solve.cpp>`__ for an example usage.

Finer-grain access
^^^^^^^^^^^^^^^^^^

Rather than simply exposing a black-box routine for applying the inverse of 
an :math:`LDL^T` or :math:`LDL^H` factorization, the following routines can 
be used to perform various solves against the lower-triangular and diagonal 
data of a frontal tree.

.. cpp:function:: void LowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )

   **TODO: More detailed description.**

.. cpp:function:: void DiagonalSolve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )

   **TODO: More detailed description.**

Data structures
---------------
The following data structures are used to help simplify the redistribution of 
right-hand side data into a form which is compatible with the distributed 
frontal tree.

DistNodalMultiVec
^^^^^^^^^^^^^^^^^

The following structure should be used when there is only a modest number of 
right-hand sides.

.. cpp:type:: struct DistNodalMultiVec<T>

   .. cpp:member:: std::vector<Matrix<T> > localNodes
   .. cpp:member:: std::vector<Matrix<T> > distNodes

   .. cpp:function:: DistNodalMultiVec( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistMultiVec<T>& X )

   .. cpp:function:: int Height() const

      Returns the length of each vector.

   .. cpp:function:: int Width() const

      Returns the number of vectors.

   .. cpp:function:: int LocalHeight() const

      Returns the total number of local rows.

.. cpp:type:: struct DistNodalMultiVec<F>

   Same as above, but this implies that the underlying datatype `F` is a field.

DistNodalMatrix
^^^^^^^^^^^^^^^

The following structure is a work in progress and should be used in cases where
there are many right-hand sides (in a sense to be made more specific later).

.. cpp:type:: struct DistNodalMatrix<T>

   .. cpp:member:: std::vector<Matrix<T> > localNodes
   .. cpp:member:: std::vector<Matrix<T> > distNodes

   .. cpp:function:: DistNodalMatrix( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistMultiVec<T>& X )

   .. cpp:function:: int Height() const

      Returns the length of each vector.

   .. cpp:function:: int Width() const

      Returns the number of vectors.

.. cpp:type:: struct DistNodalMatrix<F>

   Same as above, but this implies that the underlying datatype `F` is a field.
