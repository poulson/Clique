The DistVector class
====================
The :cpp:class:`DistVector\<T>` class is meant to...
**TODO**

.. cpp:class:: DistVector<T>

   .. rubric:: Constructors

   .. cpp:function:: DistVector()

   .. cpp:function:: DistVector( mpi::Comm comm )

   .. cpp:function:: DistVector( int height, mpi::Comm comm )

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

   .. cpp:function:: mpi::Comm Comm() const

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

      The distribution blocksize of the distributed vector:
      the process with rank ``r``'s first local row is global row
      ``r*blocksize``.

   .. cpp:function:: int FirstLocalRow() const

      The global index of the first row assigned to this process.

   .. cpp:function:: int LocalHeight() const

      The number of rows assigned to this process.

   .. rubric:: Local data

   .. cpp:function:: T GetLocal( int localRow ) const

      Returns the value of the specified local entry.

   .. cpp:function:: void SetLocal( int localRow, T value )
     
      Sets a local entry equal to a particular value.

   .. cpp:function:: void UpdateLocal( int localRow, T value )

      Add the specified value onto a local entry of the distributed vector.

   .. rubric:: For modifying the size of the vector

   .. cpp:function:: void Empty()

      Frees all resources and makes the vector have height zero.

   .. cpp:function:: void ResizeTo( int height )

      Resizes the vector to have the specified height.

   .. cpp:function:: const DistVector<T>& operator=( const DistVector<T>& x )

      Makes this vector a copy of the given vector.

.. cpp:function:: void MakeZeros( DistVector<T>& x )

.. cpp:function:: void MakeUniform( DistVector<T>& x )

.. cpp:function:: typename Base<F>::type Norm( const DistVector<F>& x )

.. cpp:function:: void Axpy( T alpha, const DistVector<T>& x, DistVector<T>& y )
