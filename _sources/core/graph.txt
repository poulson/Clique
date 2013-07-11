Graphs
======

The Graph class
---------------

The :cpp:type:`Graph` class is primarily used under-the-hood in Clique for 
graph partitioning and as the primary component of the 
:cpp:type:`SparseMatrix\<T>` class.
Nevertheless, it is possible to work directly with graphs in Clique. The usage
is extremely similar to the :cpp:type:`SparseMatrix\<T>` class, and so the 
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

.. cpp:type:: class Graph

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

   .. cpp:function:: int* SourceBuffer()
   .. cpp:function:: const int* LockedSourceBuffer() const

      Returns a (const) pointer to the contiguous array of source indices.

   .. cpp:function:: int* TargetBuffer()
   .. cpp:function:: const int* LockedTargetBuffer() const

      Returns a (const) pointer to the contiguous array of target indices.

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

The DistGraph class
-------------------
Like the :cpp:type:`Graph` class, the :cpp:type:`DistGraph` class is mainly 
used internally by Clique for (distributed) graph partitioning and for most of 
the functionality of the :cpp:type:`DistSparseMatrix\<T>` class. 
The :cpp:type:`DistGraph` class essentially provides of a subset of the 
functionality of the :cpp:type:`DistSparseMatrix\<T>` class, but an example 
of constructing the graph of a 7-point stencil over an 
:math:`n_1 \times n_2 \times n_3` grid is shown below:

.. code-block:: cpp

    DistGraph graph( N, comm );
    const int firstLocalSource = graph.FirstLocalSource();
    const int numLocalSources = graph.NumLocalSources();
    graph.StartAssembly();
    graph.Reserve( 7*numLocalSources );
    for( int iLocal=0; iLocal<numLocalSources; ++iLocal )
    {
        const int i = firstLocalSource + iLocal;
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

The first thing to notice is that all routines which relate to modifying the 
sparse matrix are wrapped with the :cpp:func:`DistGraph::StartAssembly` and 
:cpp:func:`DistGraph::StopAssembly` routines.
If the edges were all naively appended to the end of a contiguous region of 
memory, then larger and larger regions of memory would frequently need to be 
allocated and the previous contents would be copied into each new buffer. 
In order to prevent this issue, one can simply pass an upper-bound on the 
number of local updates to the :cpp:func:`DistGraph::Reserve` member function 
before inserting any edges.

After finishing assembly, all redundant edges are deleted and the edges are 
sorted lexicographically based upon their source and target indices 
(respectively).

.. cpp:type:: class DistGraph

   .. rubric:: Constructors

   .. cpp:function:: DistGraph()

      Constructs an empty graph with zero vertices over ``mpi::COMM_WORLD``.

   .. cpp:function:: DistGraph( mpi::Comm comm )

      Constructs an empty graph with zero vertices over the specified 
      communicator.

   .. cpp:function:: DistGraph( int numVertices, mpi::Comm comm )

      Constructs an empty graph with `numVertices` vertices over the specified
      communicator.

   .. cpp:function:: DistGraph( int numSources, int numTargets, mpi::Comm comm )

      Constructs an empty graph with the specified numbers of source and 
      target vertices over the given communicator.

   .. cpp:function:: DistGraph( const Graph& graph )

      Constructs a copy of the given local graph over a single-process 
      communicator.

   .. cpp:function:: DistGraph( const DistGraph& graph )

      Constructs a copy of the given distributed graph.

   .. rubric:: High-level information

   .. cpp:function:: int NumSources() const

      The number of source vertices of this graph
      (this is analogous to the height of a sparse matrix).

   .. cpp:function:: int NumTargets() const

      The number of target vertices of this graph
      (this is analogous to the width of a sparse matrix).

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

      Empties the graph and switches to the specified team of processes.

   .. cpp:function:: mpi::Comm Comm() const

      The communicator for the distributed graph.

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

      The distribution blocksize of the distributed graph: 
      the process with rank ``r``'s first local source is global source
      ``r*blocksize``.

   .. cpp:function:: int FirstLocalSource() const

      The global index of the first source assigned to this process. 

   .. cpp:function:: int NumLocalSources() const

      The number of sources assigned to this process.

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

      This should be called before inserting any edges into the graph.

   .. cpp:function:: void StopAssembly()

      This should be called after all edges have been inserted into the graph,
      as it handles converting the edges into the proper internal format.

   .. cpp:function:: void Reserve( int numLocalEdges )

      This routine should be given an upper bound on the number of edges that
      will be inserted into the graph by this process so that an appropriate
      amount of memory can be allocated to store the local edge information.

   .. cpp:function:: void Insert( int source, int target )

      Inserts an edge into our local portion of the graph which connects the 
      specified source and target vertices.

   .. cpp:function:: int Capacity() const

      The number of local edges which can be stored before a memory allocation 
      will be required (including current local edges).

   .. rubric:: Local data

   .. cpp:function:: int NumLocalEdges() const

      The number of local edges in the graph.

   .. cpp:function:: int Source( int localEdge ) const

      The source vertex of the specified local edge.

   .. cpp:function:: int Target( int localEdge ) const

      The target vertex of the specified local edge.

   .. cpp:function:: int LocalEdgeOffset( int localSource ) const

      The first local edge which begins from a local source with an 
      index greater than or equal to the given index (assuming it exists).

   .. cpp:function:: int NumConnections( int localSource ) const

      The number of edges which begin from the specified local source.

   .. cpp:function:: int* SourceBuffer()
   .. cpp:function:: const int* LockedSourceBuffer() const

      Returns a (const) pointer to the contiguous array of local source indices.

   .. cpp:function:: int* TargetBuffer()
   .. cpp:function:: const int* LockedTargetBuffer() const

      Returns a (const) pointer to the contiguous array of local target indices.

   .. rubric:: For modifying the size of the graph

   .. cpp:function:: void Empty()

      Frees all resources and modifies this graph to have zero vertices.

   .. cpp:function:: void ResizeTo( int numVertices )

      Frees all resources and modifies this graph to have ``numVertices``
      vertices.

   .. cpp:function:: void ResizeTo( int numSources, int numTargets )

      Frees all resources and modifies the graph to have the specified number
      of source and target vertices.

   .. rubric:: For copying one graph into another

   .. cpp:function:: const DistGraph& operator=( const Graph& graph )

      Sets this graph equal to the given local graph.

   .. cpp:function:: const DistGraph& operator=( const DistGraph& graph )

      Sets this graph equal to the given distributed graph.
