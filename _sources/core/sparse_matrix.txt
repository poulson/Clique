The SparseMatrix class
======================
The :cpp:class:`SparseMatrix\<T>` class is supposed to provide a simple means 
of forming sequential sparse matrices within Clique and is a simplification of 
the upcoming :cpp:class:`DistSparseMatrix\<T>` class 
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

.. cpp:class:: SparseMatrix<T>

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

   .. cpp:function:: const Graph& Graph() const

      The underlying graph of the sparse matrix.

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

   .. rubric:: For modifying the size of the matrix

   .. cpp:function:: void Empty()

      Frees all resources and sets the sparse matrix to be zero by zero.

   .. cpp:function:: void ResizeTo( int height, int width )

      Frees all resources and sets the sparse matrix to be 
      `height` :math:`\times` `width`.
