The DistSparseMatrix class
==========================
The ``DistSparseMatrix`` class is meant to provide a convenient means of 
building a distributed sparse matrix within Clique. After initializing a
real double-precision :math:`N \times N` sparse matrix distributed over the 
communicator ``comm``, e.g., with

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
sparse matrix are wrapped with ``A.StartAssembly()`` and ``A.StopAssembly()``
If the updates were all naively appended to the end of a contiguous region of 
memory, then larger and larger regions of memory would frequently need to be 
allocated and the previous contents would be copied into each new buffer.
In order to prevent this issue, one can simply pass an upper-bound on the 
number of local updates to the ``Reserve`` member function before performing 
any updates.

After finishing assembly, all updates to the same entry are combined and the 
nonzero triples, which contain the row index, column index, and value, are 
sorted in lexicographical order based upon the row and column indices 
(respectively).

.. cpp:class:: DistSparseMatrix<T>

   **TODO: Briefly describe each member variable**

   .. rubric:: Constructors

   .. cpp:function:: DistSparseMatrix()

   .. cpp:function:: DistSparseMatrix( mpi::Comm comm )

   .. cpp:function:: DistSparseMatrix( int height, mpi::Comm comm )

   .. cpp:function:: DistSparseMatrix( int height, int width, mpi::Comm comm )

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

   .. cpp:function:: int Width() const

   .. cpp:function:: const DistGraph& Graph() const

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

   .. cpp:function:: mpi::Comm Comm() const

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

   .. cpp:function:: int FirstLocalRow() const

   .. cpp:function:: int LocalHeight() const

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

   .. cpp:function:: void StopAssembly()

   .. cpp:function:: void Reserve( int numLocalEntries )

   .. cpp:function:: void Update( int row, int col, T value )

   .. cpp:function:: int Capacity() const

   .. rubric:: Local data

   .. cpp:function:: int NumLocalEntries() const

   .. cpp:function:: int Row( int localEntry ) const

   .. cpp:function:: int Col( int localEntry ) const

   .. cpp:function:: T Value( int localEntry ) const

   .. cpp:function:: int LocalEntryOffset( int localRow ) const

   .. cpp:function:: int NumConnections( int localRow ) const

   .. rubric:: For modifying the size of the matrix

   .. cpp:function:: void Empty()

   .. cpp:function:: void ResizeTo( int height, int width )

