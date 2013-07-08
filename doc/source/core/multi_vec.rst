The MultiVec class
==================
The :cpp:type:`MultiVec\<T>` class is the standard interface in Clique for
setting up multiple right-hand sides for sequential solves: it is a simplified 
version of :cpp:type:`DistMultiVec\<T>`.

.. cpp:type:: class MultiVec<T>

   .. rubric:: Constructors

   .. cpp:function:: MultiVec()

      Constructs a single vector of length zero.

   .. cpp:function:: MultiVec( int height, int width )

      Constructs `width` vectors, each of the specified height.

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

      Returns the length of the vectors.

   .. cpp:function:: int Width() const

      Returns the number of vectors.

   .. rubric:: Data

   .. cpp:function:: T Get( int row, int col ) const

      Returns the value at the specified row of a particular vector.

   .. cpp:function:: void Set( int row, int col, T value )
     
      Sets the value at the specified row of a particular vector.

   .. cpp:function:: void Update( int row, int col, T value )

      Adds a value onto the specified row entry of a particular vector.

   .. rubric:: For modifying the size of the vector

   .. cpp:function:: void Empty()

      Frees all resources and sets the multivector to a single vector of height 
      zero.

   .. cpp:function:: void ResizeTo( int height, int width )

      Reconfigures the class to have `width` vectors, each of the given height.

   .. cpp:function:: const MultiVec<T>& operator=( const MultiVec<T>& X )

      Makes this multi-vector a copy of the given multi-vector.

.. cpp:type:: class MultiVec<F>

   The same as :cpp:type:`MultiVec\<T>`, but the implication is that the 
   underlying datatype `F` is a field rather than just a ring.

.. cpp:function:: void MakeZeros( MultiVec<T>& X )

   Sets every entry in the multi-vector to zero.

.. cpp:function:: void MakeUniform( MultiVec<T>& X )

   Sets each entry in the multi-vector to a sample from the unit ball 
   appropriate for type ``T``.

.. cpp:function:: void Norms( const MultiVec<F>& X, std::vector<typename Base<F>::type>& norms )

   Returns the Euclidean norms of each vector in the multi-vector.

.. cpp:function:: typename Base<F>::type Norm( const MultiVec<F>& x )

   .. note::

      This only applies when there is only a single column.

   Returns the Euclidean norm of the column-vector ``x``. 

.. cpp:function:: void Axpy( T alpha, const MultiVec<T>& X, MultiVec<T>& Y )

   Updates :math:`Y := \alpha X + Y`.