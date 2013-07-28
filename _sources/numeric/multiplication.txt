Multiplication
==============

With sparse matrices
--------------------

.. cpp:function:: void Multiply( T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X, T beta, DistMultiVec<T>& Y )

   Forms :math:`Y := \alpha A X + \beta Y`.

With frontal trees
------------------

.. cpp:function:: void LowerMultiply( Orientation orientation, int diagOff, const DistSymmInfo& info, const DistSymmFrontTree<T>& L, DistNodalMultiVec& X )

   Overwrites :math:`X` with the lower triangle of the frontal tree multiplied
   by the original contents of :math:`X`. **TODO: More detailed description of 
   the effects of the `orientation` and `diagOff` parameters**
