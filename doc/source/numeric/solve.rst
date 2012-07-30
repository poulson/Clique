Solves
======

Solving after LDL factorization
-------------------------------

.. cpp:function:: void LDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

.. cpp:function:: void LowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

.. cpp:function:: void DiagonalSolve( const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

Solving after block LDL factorization
-------------------------------------

.. cpp:function:: void BlockLDLSolve( Orientation orientation, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

.. cpp:function:: void BlockLowerSolve( Orientation orientation, UnitOrNonUnit diag, const DistSymmInfo& info, const DistSymmFrontTree<F>& L, Matrix<F>& localX )

Data structures
---------------
**TODO**

DistNodalVector
^^^^^^^^^^^^^^^

.. cpp:type:: struct DistNodalVector<F>

   .. cpp:member:: Matrix<F> localVec

   .. cpp:function:: void Pull( const std::vector<int>& localInverseMap, const DistSymmInfo& info, const DistVector<F>& x )

   .. cpp:function:: void Push( const std::vector<int>& localInverseMap, const DistSymmInfo& info, DistVector<F>& x )

   .. cpp:function:: DistNodalVector( const std::vector<int>& localInverseMap, const DistSymmInfo& info, const DistVector<F>& x )

DistNodalMultiVector
^^^^^^^^^^^^^^^^^^^^

.. cpp:type:: struct DistNodalMultiVector<F>

   .. cpp:member:: Matrix<F> localMultiVec

   .. cpp:function:: void Pull( const std::vector<int>& localInverseMap, const DistSymmInfo& info, const DistMultiVector<F>& X )

   .. cpp:function:: void Push( const std::vector<int>& localInverseMap, const DistSymmInfo& info, DistMultiVector<F>& X )

   .. cpp:function:: DistNodalMultiVector( const std::vector<int>& localInverseMap, const DistSymmInfo& info, const DistMultiVector<F>& X )

