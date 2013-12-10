Factorization
=============

LDL factorization
-----------------

.. cpp:function:: void LDL( DistSymmInfo& info, DistSymmFrontTree<F>& L, SymmFrontType newFrontType=LDL_2D )

   Performs the specified type of symmetric or Hermitian factorization 
   (with or without intrafrontal Bunch-Kaufman pivoting, with or without 
   selective inversion, and blocked or non-blocked).
   on whether `L` is marked as Hermitian. See 
   `tests/Solve <https://github.com/poulson/Clique/blob/master/tests/Solve.cpp>`__ for an example usage.
