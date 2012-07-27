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

   **TODO: Briefly describe each member variable**

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

