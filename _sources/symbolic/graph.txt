The Graph class
===============
The ``Graph`` class is primarily used under-the-hood in Clique for graph 
partitioning and as the primary component of the ``SparseMatrix`` class.
Nevertheless, it is possible to work directly with graphs in Clique. The usage
is extremely similar to the ``SparseMatrix`` class, and so the construction
of the graph corresponding to the 7-point stencil over a 
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
``StartAssembly`` and ``StopAssembly`` member functions, to reserve space for
all of the edges in the graph with the ``Reserve`` function, and then
to add edges into the graph in various locations. 

After finishing assembly, redundant edges are deleted and the edges are sorted
lexicographically based upon the source and target indices (respectively). 
In the language of matrices, the *source* is the row index, and the *target*
is the column index.

.. cpp:class:: Graph

   .. rubric:: Constructors

   .. cpp:function:: Graph()

   .. cpp:function:: Graph( int numVertices )

   .. cpp:function:: Graph( int numSources, int numTargets )

   .. cpp:function:: Graph( const Graph& graph )

   .. cpp:function:: Graph( const DistGraph& graph )

   .. rubric:: High-level information

   .. cpp:function:: int NumSources() const

   .. cpp:function:: int NumTargets() const

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

   .. cpp:function:: void StopAssembly()

   .. cpp:function:: void Reserve( int numEdges )

   .. cpp:function:: void Insert( int source, int target )

   .. cpp:function:: int Capacity() const

   .. rubric:: Data

   .. cpp:function:: int NumEdges() const

   .. cpp:function:: int Source( int edge ) const

   .. cpp:function:: int Target( int edge ) const

   .. cpp:function:: int EdgeOffset( int source ) const

   .. cpp:function:: int NumConnections( int source ) const

   .. rubric:: For modifying the size of the graph

   .. cpp:function:: void Empty()

   .. cpp:function:: void ResizeTo( int numVertices )

   .. cpp:function:: void ResizeTo( int numSources, int numTargets )

   .. rubric:: For copying one graph into another

   .. cpp:function:: const Graph& operator=( const Graph& graph )

   .. cpp:function:: const Graph& operator=( const DistGraph& graph )

The DistGraph class
===================
Like the ``Graph`` class, the ``DistGraph`` class is mainly used internally 
by Clique for (distributed) graph partitioning and for most of the 
functionality of the ``DistSparseMatrix`` class. The ``DistGraph`` class 
essentially provides of a subset of the functionality of the 
``DistSparseMatrix`` class, but an example of constructing the graph of a 
7-point stencil over an :math:`n_1 \times n_2 \times n_3` grid is shown below:

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
sparse matrix are wrapped with the ``StartAssembly`` and ``StopAssembly`` 
routines.
If the edges were all naively appended to the end of a contiguous region of 
memory, then larger and larger regions of memory would frequently need to be 
allocated and the previous contents would be copied into each new buffer. 
In order to prevent this issue, one can simply pass an upper-bound on the 
number of local updates to the ``Reserve`` member function before inserting 
any edges.

After finishing assembly, all redundant edges are deleted and the edges are 
sorted lexicographically based upon their source and target indices 
(respectively).

.. cpp:class:: DistGraph

   .. rubric:: Constructors

   .. cpp:function:: DistGraph()

   .. cpp:function:: DistGraph( mpi::Comm comm )

   .. cpp:function:: DistGraph( int numVertices, mpi::Comm comm )

   .. cpp:function:: DistGraph( int numSources, int numTargets, mpi::Comm comm )

   .. cpp:function:: DistGraph( const Graph& graph )

   .. cpp:function:: DistGraph( const DistGraph& graph )

   .. rubric:: High-level information

   .. cpp:function:: int NumSources() const

   .. cpp:function:: int NumTargets() const

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

   .. cpp:function:: mpi::Comm Comm() const

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

   .. cpp:function:: int FirstLocalSource() const

   .. cpp:function:: int NumLocalSources() const

   .. rubric:: Assembly-related routines

   .. cpp:function:: void StartAssembly()

   .. cpp:function:: void StopAssembly()

   .. cpp:function:: void Reserve( int numLocalEdges )

   .. cpp:function:: void Insert( int source, int target )

   .. cpp:function:: int Capacity() const

   .. rubric:: Local data

   .. cpp:function:: int NumLocalEdges() const

   .. cpp:function:: int Source( int localEdge ) const

   .. cpp:function:: int Target( int localEdge ) const

   .. cpp:function:: int LocalEdgeOffset( int localSource ) const

   .. cpp:function:: int NumConnections( int localSource ) const

   .. rubric:: For modifying the size of the graph

   .. cpp:function:: void Empty()

   .. cpp:function:: void ResizeTo( int numVertices )

   .. cpp:function:: void ResizeTo( int numSources, int numTargets )

   .. rubric:: For copying one graph into another

   .. cpp:function:: const DistGraph& operator=( const Graph& graph )

   .. cpp:function:: const DistGraph& operator=( const DistGraph& graph )

