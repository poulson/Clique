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
