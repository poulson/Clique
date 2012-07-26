The SparseMatrix class
======================
The ``SparseMatrix`` class is supposed to provide a simple means of forming
sequential sparse matrices within Clique and is a simplification of the 
upcoming ``DistSparseMatrix`` class (which most users will work with).
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
``StartAssembly`` and ``StopAssembly`` member functions, to reserve space for
all of the updates to the sparse matrix with the ``Reserve`` function, and then
to add entries onto the sparse matrix in various locations. 

After finishing assembly, all updates to the same entry are combined and the 
nonzero triples, which contain the row index, column index, and value, are 
sorted in lexicographical order based upon the row and column indices 
(respectively).

.. cpp:class:: SparseMatrix<T>

   **TODO: Briefly describe each member variable**

   .. rubric:: Constructors

   .. cpp:function:: SparseMatrix()

   .. cpp:function:: SparseMatrix( int height )

   .. cpp:function:: SparseMatrix( int height, int width )

   .. rubric:: High-level information

   .. cpp:function:: int Height() const

   .. cpp:function:: int Width() const

   .. cpp:function:: const Graph& Graph() const

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

   .. cpp:function:: void StopAssembly()

   .. cpp:function:: void Reserve( int numEntries )

   .. cpp:function:: void Update( int row, int col, T value )

   .. cpp:function:: int Capacity() const

   .. rubric:: Data

   .. cpp:function:: int NumEntries() const

   .. cpp:function:: int Row( int entry ) const

   .. cpp:function:: int Col( int entry ) const

   .. cpp:function:: T Value( int entry ) const

   .. cpp:function:: int EntryOffset( int row ) const

   .. cpp:function:: int NumConnections( int row ) const

   .. rubric:: For modifying the size of the matrix

   .. cpp:function:: void Empty()

   .. cpp:function:: void ResizeTo( int height, int width )

