Solves
======

Solving after LDL factorization
-------------------------------

.. cpp:function:: void LDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

   After having performed an :math:`LDL^T` or :math:`LDL^H` factorization 
   (e.g., via :cpp:func:`LDL`), this routine can be used to solve a set of 
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

Solving after block LDL factorization
-------------------------------------

.. cpp:function:: void BlockLDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

   After having performed a block :math:`LDL^T` or :math:`LDL^H` factorization 
   (e.g., via :cpp:func:`BlockLDL`), this routine can be used to solve a set of 
   right-hand sides. Note that `localX` can be easily generated using the 
   :cpp:type:`DistNodalMultiVector\<F>` class. See
   `tests/MultiVectorSolve <https://github.com/poulson/Clique/blob/master/tests/MultiVectorSolve.cpp>`__ for an example usage with the :cpp:func:`LDL` routine.

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
