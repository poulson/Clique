The Graph class
===============
The :cpp:class:`Graph` class is primarily used under-the-hood in Clique for 
graph partitioning and as the primary component of the 
:cpp:class:`SparseMatrix\<T>` class.
Nevertheless, it is possible to work directly with graphs in Clique. The usage
is extremely similar to the :cpp:class:`SparseMatrix\<T>` class, and so the 
construction of the graph corresponding to the 7-point stencil over a 
:math:`n_1 \times n_2 \times n_3` grid is briefly demonstrated below:

.. code-block:: cpp

   Graph graph( N );
   graph.StartAssembly();
   graph.Reserve( 7*N );
   for( int i=0; i<N; ++i )
   {
       const int x = i % n1;
       const int y = (i/n1) % n2;
       const int z = i/(n1*n2);
       
       graph.Insert( i, i );
       if( x != 0 )
           graph.Insert( i, i-1 );
       if( x != n1-1 )
           graph.Insert( i, i+1 );
       if( y != 0 )
           graph.Insert( i, i-n1 );
       if( y != n2-1 )
           graph.Insert( i, i+n1 );
       if( z != 0 )
           graph.Insert( i, i-n1*n2 );
       if( z != n3-1 )
           graph.Insert( i, i+n1*n2 );
   }
   graph.StopAssembly();

The basic pattern is to wrap the entire assembly process with the 
:cpp:func:`Graph::StartAssembly` and :cpp:func:`Graph::StopAssembly` member 
functions, to reserve space for
all of the edges in the graph with the :cpp:func:`Graph::Reserve` function, 
and then to add edges into the graph in various locations. 

After finishing assembly, redundant edges are deleted and the edges are sorted
lexicographically based upon the source and target indices (respectively). 
In the language of matrices, the *source* is the row index, and the *target*
is the column index.

.. cpp:class:: Graph

   .. rubric:: Constructors

   .. cpp:function:: Graph()

      Constructs an empty graph with zero vertices.

   .. cpp:function:: Graph( int numVertices )

      Constructs an empty graph with `numVertices` vertices.

   .. cpp:function:: Graph( int numSources, int numTargets )

      Constructs an empty graph with `numSources` source vertices and 
      `numTargets` target vertices: this is analogous to creating an empty
      `numSources` :math:`\times` `numTargets` sparse matrix.

   .. cpp:function:: Graph( const Graph& graph )

      Constructs a copy of the given graph.

   .. cpp:function:: Graph( const DistGraph& graph )

      .. note::
         
         The distributed graph can only be distributed over a single process.

      Converts a trivial distributed graph into an explicitly local graph.

   .. rubric:: High-level information

   .. cpp:function:: int NumSources() const

      The number of source vertices in the graph; this is analogous to the
      height of a sparse matrix.

   .. cpp:function:: int NumTargets() const

      The number of target vertices in the graph; this is analogous to the
      width of a sparse matrix.

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

      This should be called before inserting any edges into the graph.

   .. cpp:function:: void StopAssembly()

      This should be called after all edges have been inserted into the graph,
      as it handles converting the edges into the proper internal format.

   .. cpp:function:: void Reserve( int numEdges )

      This routine should be given an upper bound on the number of edges that
      will be inserted into the graph so that a sufficient amount of memory 
      can be allocated to store all of the edge information.

   .. cpp:function:: void Insert( int source, int target )

      Inserts an edge starting from the specified source vertex which ends
      at the target vertex.

   .. cpp:function:: int Capacity() const

      The total number of edges which can be stored before a memory allocation
      (including current edges).

   .. rubric:: Data

   .. cpp:function:: int NumEdges() const

      The number of edges which have been inserted into the graph.

   .. cpp:function:: int Source( int edge ) const

      The source vertex of the specified edge.

   .. cpp:function:: int Target( int edge ) const

      The target vertex of the specified edge.

   .. cpp:function:: int EdgeOffset( int source ) const

      The first edge which begins at a source vertex with an equal or greater
      index (assuming it exists).

   .. cpp:function:: int NumConnections( int source ) const

      The number of edges which start from the specified source vertex.

   .. rubric:: For modifying the size of the graph

   .. cpp:function:: void Empty()

      Frees all resources and modified this graph to have zero vertices.

   .. cpp:function:: void ResizeTo( int numVertices )

      Frees all resources and modifies this graph to have `numVertices` 
      vertices.

   .. cpp:function:: void ResizeTo( int numSources, int numTargets )

      Frees all resources and modifies this graph to have `numSources`
      source vertices and `numTargets` target vertices.

   .. rubric:: For copying one graph into another

   .. cpp:function:: const Graph& operator=( const Graph& graph )

      Sets this graph equal to the given graph.

   .. cpp:function:: const Graph& operator=( const DistGraph& graph )

      .. note::

         The distributed graph can only be distributed over a single process.

      Sets this graph equal to the given graph.
