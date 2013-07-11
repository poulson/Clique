Sparse matrices
===============

The SparseMatrix class
----------------------
The :cpp:type:`SparseMatrix\<T>` class is supposed to provide a simple means 
of forming sequential sparse matrices within Clique and is a simplification of 
the upcoming :cpp:type:`DistSparseMatrix\<T>` class 
(which most users will work with).
One might start by first constructing an empty double-precision 
:math:`N \times N` sparse matrix, e.g., with

.. code-block:: cpp

    SparseMatrix<double> A( N );

A typical assembly of a sparse matrix (i.e., using a 7-point finite-difference stencil over an :math:`n_1 \times n_2 \times n_3` grid) might be

.. code-block:: cpp

   A.StartAssembly();
   A.Reserve( 7*N );
   for( int i=0; i<N; ++i )
   {
       const int x = i % n1;
       const int y = (i/n1) % n2;
       const int z = i/(n1*n2);
       
       A.Update( i, i, 6. );
       if( x != 0 )
           A.Update( i, i-1, -1. );
       if( x != n1-1 )
           A.Update( i, i+1, -1. );
       if( y != 0 )
           A.Update( i, i-n1, -1. );
       if( y != n2-1 )
           A.Update( i, i+n1, -1. );
       if( z != 0 )
           A.Update( i, i-n1*n2, -1. );
       if( z != n3-1 )
           A.Update( i, i+n1*n2, -1. );
   }
   A.StopAssembly();

The basic pattern is to wrap the entire assembly process with the 
:cpp:func:`SparseMatrix\<T>::StartAssembly` and 
:cpp:func:`SparseMatrix\<T>::StopAssembly` member functions, 
to reserve space for all of the updates to the sparse matrix with the 
:cpp:func:`SparseMatrix\<T>::Reserve` function, and then
to add entries onto the sparse matrix in various locations. 

After finishing assembly, all updates to the same entry are combined and the 
nonzero triples, which contain the row index, column index, and value, are 
sorted in lexicographical order based upon the row and column indices 
(respectively).

.. cpp:type:: class SparseMatrix<T>

   .. rubric:: Constructors

   .. cpp:function:: SparseMatrix()

      Constructs an empty :math:`0 \times 0` sparse matrix.

   .. cpp:function:: SparseMatrix( int height )

      Constructs an empty `height` :math:`\times` `height` sparse matrix.

   .. cpp:function:: SparseMatrix( int height, int width )

      Constructs an empty `height` :math:`\times` `width` sparse matrix.

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

      The height of the sparse matrix.

   .. cpp:function:: int Width() const

      The width of the sparse matrix.

   .. cpp:function:: Graph& Graph()
   .. cpp:function:: const Graph& LockedGraph() const

      The underlying (immutable) graph of the sparse matrix.

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

      This should be called before updating any entries of the sparse matrix.

   .. cpp:function:: void StopAssembly()

      This should be called after all updates have been applied to the sparse
      matrix, as it handles combining updates of the same entry and then
      sorting entry data into the proper internal format.

   .. cpp:function:: void Reserve( int numEntries )

      This routine should be given an upper bound on the number of updates 
      that will be applied to the sparse matrix so that a sufficient amount 
      of memory can be allocated to store all of the update information.

   .. cpp:function:: void Update( int row, int col, T value )

      Adds the specified value to the entry of the sparse matrix with the 
      given row and column indices.

   .. cpp:function:: int Capacity() const

      The number of updates that can be applied to the sparse matrix before
      a memory allocation (including current updates).

   .. rubric:: Data

   .. cpp:function:: int NumEntries() const

      The number of nonzero entries in the sparse matrix.

   .. cpp:function:: int Row( int entry ) const

      The row index of the given nonzero entry.

   .. cpp:function:: int Col( int entry ) const

      The column index of the given nonzero entry.

   .. cpp:function:: T Value( int entry ) const

      The numerical value of the given nonzero entry.

   .. cpp:function:: int EntryOffset( int row ) const

      The first nonzero entry within a row with index greater than or equal
      to the given value (assuming it exists).

   .. cpp:function:: int NumConnections( int row ) const

      The number of nonzero entries in the specified row.

   .. cpp:function:: int* SourceBuffer()
   .. cpp:function:: const int* LockedSourceBuffer() const

      Returns a (const) pointer to the contiguous array of source indices.

   .. cpp:function:: int* TargetBuffer()
   .. cpp:function:: const int* LockedTargetBuffer() const

      Returns a (const) pointer to the contiguous array of target indices.

   .. cpp:function:: T* ValueBuffer()
   .. cpp:function:: const T* LockedValueBuffer() const

      Returns a (const) pointer to the contiguous array of nonzero values.

   .. rubric:: For modifying the size of the matrix

   .. cpp:function:: void Empty()

      Frees all resources and sets the sparse matrix to be zero by zero.

   .. cpp:function:: void ResizeTo( int height, int width )

      Frees all resources and sets the sparse matrix to be 
      `height` :math:`\times` `width`.

.. cpp:type:: class SparseMatrix<F>

   This is the same as :cpp:type:`SparseMatrix\<T>`, but the implication is that
   the underlying datatype `F` is a field rather than just a ring.

The DistSparseMatrix class
--------------------------
The :cpp:type:`DistSparseMatrix\<T>` class is meant to provide a convenient 
means of building a distributed sparse matrix within Clique. 
After initializing a real double-precision :math:`N \times N` sparse matrix 
distributed over the communicator ``comm``, e.g., with

.. code-block:: cpp

    DistSparseMatrix<double> A( N, comm );

the first row that a particular process owns, as well as how many rows it
owns, can be returned with 

.. code-block:: cpp

    const int firstLocalRow = A.FirstLocalRow();
    const int localHeight = A.LocalHeight();

Each process is then responsible for specifying the nonzeros in rows 
``[firstLocalRow,firstLocalRow+localHeight)``, and this is accomplished by 
*updating* particular locations in the matrix. It is instructive to look at an
example of filling a 7-point finite-difference stencil over an :math:`n_1 \times n_2 \times n_3` grid before describing the individual functions

.. code-block:: cpp

    A.StartAssembly();
    A.Reserve( 7*localHeight );
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const int i = firstLocalRow + iLocal;
        const int x = i % n1;
        const int y = (i/n1) % n2;
        const int z = i/(n1*n2);

        A.Update( i, i, 6. );
        if( x != 0 )
            A.Update( i, i-1, -1. );
        if( x != n1-1 )
            A.Update( i, i+1, -1. );
        if( y != 0 )
            A.Update( i, i-n1, -1. );
        if( y != n2-1 )
            A.Update( i, i+n1, -1. );
        if( z != 0 )
            A.Update( i, i-n1*n2, -1. );
        if( z != n3-1 )
            A.Update( i, i+n1*n2, -1. );
    }
    A.StopAssembly();

The first thing to notice is that all routines which relate to modifying the 
sparse matrix are wrapped with :cpp:func:`DistSparseMatrix\<T>::StartAssembly` 
and :cpp:func:`DistSparseMatrix\<T>::StopAssembly`.
If the updates were all naively appended to the end of a contiguous region of 
memory, then larger and larger regions of memory would frequently need to be 
allocated and the previous contents would be copied into each new buffer.
In order to prevent this issue, one can simply pass an upper-bound on the 
number of local updates to the :cpp:func:`DistSparseMatrix\<T>::Reserve` 
member function before performing any updates.

After finishing assembly, all updates to the same entry are combined and the 
nonzero triples, which contain the row index, column index, and value, are 
sorted in lexicographical order based upon the row and column indices 
(respectively).

.. cpp:type:: class DistSparseMatrix<T>

   .. rubric:: Constructors

   .. cpp:function:: DistSparseMatrix()

      Constructs an empty :math:`0 \times 0` sparse matrix over 
      ``mpi::COMM_WORLD``.

   .. cpp:function:: DistSparseMatrix( mpi::Comm comm )

      Constructs an empty :math:`0 \times 0` sparse matrix over the specified
      communicator.

   .. cpp:function:: DistSparseMatrix( int height, mpi::Comm comm )

      Constructs a `height` :math:`\times` `height` sparse matrix over the 
      specified communicator.

   .. cpp:function:: DistSparseMatrix( int height, int width, mpi::Comm comm )

      Constructs a `height` :math:`\times` `width` sparse matrix over the 
      specified communicator.

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

      The height of the sparse matrix.

   .. cpp:function:: int Width() const

      The width of the sparse matrix.

   .. cpp:function:: DistGraph& DistGraph()
   .. cpp:function:: const DistGraph& LockedDistGraph() const

      The underlying (immutable) distributed graph of the sparse matrix.

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

      Empties the sparse matrix and switches to the specified team of 
      processes.

   .. cpp:function:: mpi::Comm Comm() const

      The communicator for the distributed sparse matrix. 

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

      The distribution blocksize of the distributed sparse matrix:
      the process with rank ``r``'s first local row is global row
      ``r*blocksize``.

   .. cpp:function:: int FirstLocalRow() const

      The global index of the first row assigned to this process.

   .. cpp:function:: int LocalHeight() const

      The number of rows assigned to this process.

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

      This should be called before applying any updates to the sparse matrix.

   .. cpp:function:: void StopAssembly()

      This should be called after all updates have been applied to the sparse
      matrix, as it handles combining updates to the same entry and converting
      the entry information into the proper internal format.

   .. cpp:function:: void Reserve( int numLocalEntries )

      This routine should be given an upper bound on the number of updates 
      that will be applied by this process so that an appropriate amount of 
      memory can be allocated to store the update information.

   .. cpp:function:: void Update( int row, int col, T value )

      Add the specified value onto the entry of our local portion of the 
      sparse matrix with the specified indices. 

   .. cpp:function:: int Capacity() const

      The number of updates which can be applied before a memory allocation
      will be required (including current local updates).

   .. rubric:: Local data

   .. cpp:function:: int NumLocalEntries() const

      The number of local nonzero entries in the sparse matrix.

   .. cpp:function:: int Row( int localEntry ) const

      The row index of the specified local nonzero entry.

   .. cpp:function:: int Col( int localEntry ) const

      The column index of the specified local nonzero entry.

   .. cpp:function:: T Value( int localEntry ) const

      The numerical value of the specified local nonzero entry.

   .. cpp:function:: int LocalEntryOffset( int localRow ) const

      The first local nonzero entry which lies within a row with local index 
      greater than or equal to the specified value (assuming it exists).

   .. cpp:function:: int NumConnections( int localRow ) const

      The number of nonzeros within the specified local row.

   .. cpp:function:: int* SourceBuffer()
   .. cpp:function:: const int* LockedSourceBuffer() const

      Returns a (const) pointer to the contiguous array of local source indices.

   .. cpp:function:: int* TargetBuffer()
   .. cpp:function:: const int* LockedTargetBuffer() const

      Returns a (const) pointer to the contiguous array of local target indices.

   .. cpp:function:: T* ValueBuffer()
   .. cpp:function:: const T* LockedValueBuffer() const

      Returns a (const) pointer to the contiguous array of local nonzero values.

   .. rubric:: For modifying the size of the matrix

   .. cpp:function:: void Empty()

      Frees all resources and sets the sparse matrix to be zero by zero.

   .. cpp:function:: void ResizeTo( int height, int width )

      Frees all resources and sets the sparse matrix to be 
      `height` :math:`\times` `width`.

.. cpp:type:: class DistSparseMatrix<F>

   This is the same as :cpp:type:`DistSparseMatrix\<T>`, but the implication 
   is that the underlying datatype `F` is a field rather than just a ring.
