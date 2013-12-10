Solves
======

Black-box solvers
-----------------

.. cpp:function:: void SymmetricSolve( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool conjugate=false, bool sequential=true, int numDistSeps=1, int numSeqSeps=1, int cutoff=128 )

   Overwrites :math:`x` with :math:`A^{-1} x`, where :math:`A` is assumed to be 
   symmetric (or Hermitian if `conjugate` is `true`). 
   The main optional argument is 
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
   Hermitian. This is simply a wrapper to `SymmetricSolve` with `conjugate`
   set to `true`.

Solving after factorization
---------------------------

.. cpp:function:: void Solve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )

   After factorization (via :cpp:func:`LDL`), this routine can be used to 
   solve a set of right-hand sides. See
   `tests/Solve <https://github.com/poulson/Clique/blob/master/tests/Solve.cpp>`__ for an example usage.

Finer-grain access
^^^^^^^^^^^^^^^^^^

Rather than simply exposing a black-box routine for applying the inverse of 
a factorization, the following routines can be used to perform various solves 
against the lower-triangular and (quasi-)diagonal data of a frontal tree.

.. cpp:function:: void LowerSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )

   **TODO: More detailed description.**

.. cpp:function:: void DiagonalSolve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, DistNodalMultiVec<F>& X )

   **TODO: More detailed description.**
