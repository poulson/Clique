Sparsity-independent vectors
============================

The MultiVec class
------------------
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

The DistMultiVec class
----------------------
The :cpp:type:`DistMultiVec\<T>` class is the standard interface in Clique 
for setting up several right-hand sides for solves: it uses the same simple 
one-dimensional distribution as :cpp:type:`DistSparseMatrix\<T>` and is meant 
to be simple to use.

.. cpp:type:: class DistMultiVec<T>

   .. rubric:: Constructors

   .. cpp:function:: DistMultiVec()

      Constructs a single vector of length zero over ``mpi::COMM_WORLD``.

   .. cpp:function:: DistMultiVec( mpi::Comm comm )

      Constructs a single vector of length zero over the specified communicator.

   .. cpp:function:: DistMultiVec( int height, int width, mpi::Comm comm )

      Constructs `width` vectors, each of the given height, over a particular 
      communicator.

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

      The length of the vectors.

   .. cpp:function:: int Width() const

      The number of vectors.

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

      Empties the multi-vector and reconfigures it to be distributed over the 
      specified communicator.

   .. cpp:function:: mpi::Comm Comm() const

      The communciator for the distributed multi-vector.

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

      The distribution blocksize of the multi-vector: the process with 
      rank ``r``'s first local row is global row ``r*blocksize``.

   .. cpp:function:: int FirstLocalRow() const

      The global index of the first row assigned to this process.

   .. cpp:function:: int LocalHeight() const

      The number of rows assigned to this process.

   .. rubric:: Local data

   .. cpp:function:: T GetLocal( int localRow, int col ) const

      The value of the specified local entry in the `col` vector.

   .. cpp:function:: void SetLocal( int localRow, int col, T value )
     
      Sets a local entry equal to a particular value.

   .. cpp:function:: void UpdateLocal( int localRow, int col, T value )

      Add the specified value onto a local entry of the distributed 
      multi-vector.

   .. rubric:: For modifying the size of the multi-vector

   .. cpp:function:: void Empty()

      Frees all resources and empties the multi-vector.

   .. cpp:function:: void ResizeTo( int height, int width )

      Resizes the multi-vector to have `width` vectors, each of the specified
      height.

   .. cpp:function:: const DistMultiVec<T>& operator=( const DistMultiVec<T>& X )

      Makes this multi-vector a copy of the given multi-vector.

.. cpp:type:: class DistMultiVec<F>

   The same as :cpp:type:`DistMultiVec\<T>`, but the implication is that the
   underlying datatype `F` is a field rather than just a ring.

.. cpp:function:: void MakeZeros( DistMultiVec<T>& X )

   Sets every entry in the multi-vector to zero.

.. cpp:function:: void MakeUniform( DistMultiVec<T>& X )

   Sets each entry in the multi-vector to a sample from the unit ball 
   appropriate for type ``T``.

.. cpp:function:: void Norms( const DistMultiVec<F>& X, std::vector<typename Base<F>::type>& norms )

   Returns the Euclidean norms of each of the column vectors.

.. cpp:function:: typename Base<F>::type Norm( const DistMultiVec<F>& x )

   .. note::

      This only applies when there is only a single column.
   
   Returns the Euclidean norm of the column-vector ``x``. 

.. cpp:function:: void Axpy( T alpha, const DistMultiVec<T>& X, DistMultiVec<T>& Y )

   Updates :math:`Y := \alpha X + Y`.
