The DistVector class
====================
The :cpp:type:`DistVector\<T>` class is the standard interface in Clique for
setting up right-hand sides for solves: it uses the same simple one-dimensional
distribution as :cpp:type:`DistSparseMatrix\<T>` and is meant to be simple to
use.

.. cpp:type:: class DistVector<T>

   .. rubric:: Constructors

   .. cpp:function:: DistVector()

      Constructs a vector of length zero over ``mpi::COMM_WORLD``.

   .. cpp:function:: DistVector( mpi::Comm comm )

      Constructs a vector of length zero over the specified communicator.

   .. cpp:function:: DistVector( int height, mpi::Comm comm )

      Constructs a vector of the given height over a particular communicator.

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

      The length of the vector.

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

      Reconfigures the vector to be of length zero and distributed over the 
      specified communicator.

   .. cpp:function:: mpi::Comm Comm() const

      The communciator for the distributed vector.

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

      The distribution blocksize of the vector: the process with rank ``r``'s
      first local row is global row ``r*blocksize``.

   .. cpp:function:: int FirstLocalRow() const

      The global index of the first row assigned to this process.

   .. cpp:function:: int LocalHeight() const

      The number of rows assigned to this process.

   .. rubric:: Local data

   .. cpp:function:: T GetLocal( int localRow ) const

      The value of the specified local entry.

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

.. cpp:type:: class DistVector<F>

   The same as :cpp:type:`DistVector\<T>`, but the implication is that the 
   underlying datatype `F` is a field rather than just a ring.

.. cpp:function:: void MakeZeros( DistVector<T>& x )

   Sets every entry in the vector to zero.

.. cpp:function:: void MakeUniform( DistVector<T>& x )

   Sets each entry in the vector to a sample from the unit ball appropriate 
   for type ``T``.

.. cpp:function:: typename Base<F>::type Norm( const DistVector<F>& x )

   The Euclidean norm of the vector.

.. cpp:function:: void Axpy( T alpha, const DistVector<T>& x, DistVector<T>& y )

   Updates :math:`y := \alpha x + y`.
