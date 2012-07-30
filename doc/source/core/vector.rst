The Vector class
================
The :cpp:type:`Vector\<T>` class is the standard interface in Clique for
setting up right-hand sides for sequential solves: it is a simplified version of
:cpp:type:`DistVector\<T>`.

.. cpp:type:: class Vector<T>

   .. rubric:: Constructors

   .. cpp:function:: Vector()

      Constructs a vector of length zero.

   .. cpp:function:: Vector( int height )

      Constructs a vector of the specified height.

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

      The length of the vector.

   .. rubric:: Data

   .. cpp:function:: T Get( int row ) const

      The value of the specified entry.

   .. cpp:function:: void Set( int row, T value )
     
      Sets an entry equal to a particular value.

   .. cpp:function:: void Update( int row, T value )

      Add the specified value onto an entry.

   .. rubric:: For modifying the size of the vector

   .. cpp:function:: void Empty()

      Frees all resources and makes the vector have height zero.

   .. cpp:function:: void ResizeTo( int height )

      Resizes the vector to have the specified height.

   .. cpp:function:: const Vector<T>& operator=( const Vector<T>& x )

      Makes this vector a copy of the given vector.

.. cpp:type:: class Vector<F>

   The same as :cpp:type:`Vector\<T>`, but the implication is that the 
   underlying datatype `F` is a field rather than just a ring.

.. cpp:function:: void MakeZeros( Vector<T>& x )

   Sets every entry in the vector to zero.

.. cpp:function:: void MakeUniform( Vector<T>& x )

   Sets each entry in the vector to a sample from the unit ball appropriate 
   for type ``T``.

.. cpp:function:: typename Base<F>::type Norm( const Vector<F>& x )

   The Euclidean norm of the vector.

.. cpp:function:: void Axpy( T alpha, const Vector<T>& x, Vector<T>& y )

   Updates :math:`y := \alpha x + y`.
