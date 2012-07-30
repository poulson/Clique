Multiplication
==============

With sparse matrices
--------------------

.. cpp:function:: void Multiply( T alpha, const DistSparseMatrix<T>& A, const DistVector<T>& x, T beta, DistVector<T>& y )

   Forms :math:`y := \alpha A x + \beta y`.

.. cpp:function:: void Multiply( T alpha, const DistSparseMatrix<T>& A, const DistMultiVector<T>& X, T beta, DistMultiVector<T>& Y )

   Forms :math:`Y := \alpha A X + \beta Y`.

With frontal trees
------------------

.. cpp:function:: void LowerMultiply( Orientation orientation, UnitOrNonUnit diag, int diagOffset, const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )

   Overwrites :math:`X` with the lower triangle of the frontal tree multiplied
   by the original contents of :math:`X`. **TODO: More detailed description of 
   the effects of the `orientation`, `diag`, and `diagOffset` parameters**
