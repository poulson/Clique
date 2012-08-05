Solves
======

Black-box solvers
-----------------

.. cpp:function:: void SymmetricSolve( const DistSparseMatrix<F>& A, DistVector<F>& x, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )
.. cpp:function:: void SymmetricSolve( const DistSparseMatrix<F>& A, DistMultiVector<F>& X, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )

   Overwrites :math:`x` with :math:`A^{-1} x`, where :math:`A` is assumed to be 
   symmetric. The main optional argument is 
   whether or not graph partitioning should be performed sequentially or in 
   parallel. For modest-to-large numbers of processes, one should only use 
   parallel graph partitioning if the entire graph for the sparse matrix cannot
   be stored on a single process. The remaining optional arguments determine 
   how many distributed and sequential separators should be tested for each 
   separator (`numDistSeps` and `numSeqSeps`) and what the maximum allowed 
   subdomain size is (`cutoff`). See
   `tests/SimpleVectorSolve <https://github.com/poulson/Clique/blob/master/tests/SimpleVectorSolve.cpp>`__ and 
   `tests/SimpleMultiVectorSolve <https://github.com/poulson/Clique/blob/master/tests/SimpleMultiVectorSolve.cpp>`__ for an example usages.

.. cpp:function:: void HermitianSolve( const DistSparseMatrix<F>& A, DistVector<F>& x, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )
.. cpp:function:: void HermitianSolve( const DistSparseMatrix<F>& A, DistMultiVector<F>& X, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )

   Overwrites :math:`x` with :math:`A^{-1} x`, where :math:`A` is assumed to be
   Hermitian. The optional arguments are identical to those of 
   :cpp:func:`SymmetricSolve`.

Solving after factorization
---------------------------

.. cpp:function:: void Solve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

   After having performed an :math:`LDL^T` or :math:`LDL^H` factorization 
   (via :cpp:func:`LDL`), this routine can be used to solve a set of 
   right-hand sides. Note that `localX` can be easily generated using the 
   :cpp:type:`DistNodalMultiVector\<F>` class. See
   `tests/MultiVectorSolve <https://github.com/poulson/Clique/blob/master/tests/MultiVectorSolve.cpp>`__ for an example usage.

Finer-grain access
^^^^^^^^^^^^^^^^^^

Rather than simply exposing a black-box routine for applying the inverse of 
an :math:`LDL^T` or :math:`LDL^H` factorization, the following routines can 
be used to perform various solves against the lower-triangular and diagonal 
data of a frontal tree.

.. cpp:function:: void LowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

   **TODO: More detailed description.**

.. cpp:function:: void DiagonalSolve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

   **TODO: More detailed description.**

Finer-grain access
^^^^^^^^^^^^^^^^^^

The following routine is an analogue to :cpp:func:`LowerSolve` for the case 
where a block :math:`LDL^T` or :math:`LDL^H` factorization was formed instead 
of a standard one.

.. cpp:function:: void BlockLowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

Data structures
---------------
The following data structures are used to help simplify the redistribution of 
right-hand side data into a form which is compatible with the distributed 
frontal tree.

DistNodalVector
^^^^^^^^^^^^^^^

.. cpp:type:: struct DistNodalVector<T>

   .. cpp:member:: Matrix<T> localVec

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistVector<T>& x )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistVector<T>& x )

   .. cpp:function:: DistNodalVector( const DistMap& inverseMap, const DistSymmInfo& info, const DistVector<T>& x )

.. cpp:type:: struct DistNodalVector<F>

   Same as above, but this implies that the underlying datatype `F` is a field.

DistNodalMultiVector
^^^^^^^^^^^^^^^^^^^^

.. cpp:type:: struct DistNodalMultiVector<T>

   .. cpp:member:: Matrix<T> localMultiVec

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVector<T>& X )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistMultiVector<T>& X )

   .. cpp:function:: DistNodalMultiVector( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVector<T>& X )

.. cpp:type:: struct DistNodalMultiVector<F>

   Same as above, but this implies that the underlying datatype `F` is a field.
