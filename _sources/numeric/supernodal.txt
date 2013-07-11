Sparsity-dependent vectors
==========================
The following data structures are used to help simplify the redistribution of 
right-hand side data into a form which is compatible with the distributed 
frontal tree.

DistNodalMultiVec
-----------------

The following structure should be used when there is only a modest number of 
right-hand sides.

.. cpp:type:: struct DistNodalMultiVec<T>

   .. cpp:member:: std::vector<Matrix<T> > localNodes
   .. cpp:member:: std::vector<DistMatrix<T,VC,STAR> > distNodes

   .. cpp:function:: DistNodalMultiVec( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: DistNodalMultiVec( const DistNodalMatrix<T>& X )

   .. cpp:function:: const DistNodalMultiVec<T>& operator=( const DistNodalMatrix<T>& X )

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistMultiVec<T>& X )

   .. cpp:function:: int Height() const

      Returns the length of each vector.

   .. cpp:function:: int Width() const

      Returns the number of vectors.

   .. cpp:function:: int LocalHeight() const

      Returns the total number of local rows.

.. cpp:type:: struct DistNodalMultiVec<F>

   Same as above, but this implies that the underlying datatype `F` is a field.

DistNodalMatrix
---------------

The following structure should be used when there are many right-hand sides.

.. cpp:type:: struct DistNodalMatrix<T>

   .. cpp:member:: std::vector<Matrix<T> > localNodes
   .. cpp:member:: std::vector<DistMatrix<T> > distNodes

   .. cpp:function:: DistNodalMatrix( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: DistNodalMatrix( const DistNodalMultiVec<T>& X )

   .. cpp:function:: const DistNodalMatrix<T>& operator=( const DistNodalMultiVec<T>& X )

   .. cpp:function:: void Pull( const DistMap& inverseMap, const DistSymmInfo& info, const DistMultiVec<T>& X )

   .. cpp:function:: void Push( const DistMap& inverseMap, const DistSymmInfo& info, DistMultiVec<T>& X )

   .. cpp:function:: int Height() const

      Returns the length of each vector.

   .. cpp:function:: int Width() const

      Returns the number of vectors.

   .. cpp:member:: mutable std::vector<MatrixCommMeta> commMetas
   .. cpp:function:: void ComputeCommMetas( const DistSymmInfo& info ) const

.. cpp:type:: struct DistNodalMatrix<F>

   Same as above, but this implies that the underlying datatype `F` is a field.
