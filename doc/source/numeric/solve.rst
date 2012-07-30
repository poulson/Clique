Solves
======

Solving after LDL factorization
-------------------------------

.. cpp:function:: void LDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )

.. cpp:function:: void LowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )

.. cpp:function:: void DiagonalSolve( const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )

Solving after block LDL factorization
-------------------------------------

.. cpp:function:: void BlockLDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )

.. cpp:function:: void BlockLowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<T>& L, Matrix<T>& localX )

Data structures
---------------
**TODO**

DistNodalVector
^^^^^^^^^^^^^^^

.. cpp:type:: struct DistNodalVector<T>

   .. cpp:member:: Matrix<T> localVec

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistVector<T>& x )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistVector<T>& x )

   .. cpp:function:: DistNodalVector( const DistMap& inverseMap, const DistSymmInfo& info, const DistVector<T>& x )

DistNodalMultiVector
^^^^^^^^^^^^^^^^^^^^

.. cpp:type:: struct DistNodalMultiVector<T>

   .. cpp:member:: Matrix<T> localMultiVec

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVector<T>& X )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistMultiVector<T>& X )

   .. cpp:function:: DistNodalMultiVector( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVector<T>& X )

