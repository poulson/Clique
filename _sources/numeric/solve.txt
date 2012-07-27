Solving after factorization
===========================

After LDL factorization
-----------------------

.. cpp:function:: void LDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

.. cpp:function:: void LowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

.. cpp:function:: void DiagonalSolve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

After block LDL factorization
-----------------------------

.. cpp:function:: void BlockLDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

.. cpp:function:: void BlockLowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

