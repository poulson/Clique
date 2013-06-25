Input/output
************

Display
-------

.. cpp:function:: void Display( const Graph& graph, std::string title="Graph" )
.. cpp:function:: void Display( const DistGraph& graph, std::string title="DistGraph" )

   Visualizes the edges in the graph 
   (if Qt5 is not available, ``Print`` instead).

.. cpp:function:: void Display( const SparseMatrix<T>& A, std::string title="SparseMatrix" )
.. cpp:function:: void Display( const DistSparseMatrix<T>& A, std::string title="DistSparseMatrix" )

   Visualize the entries in the sparse matrix
   (if Qt5 is not available, ``Print`` instead).

DisplayLocal
------------

.. cpp:function:: void DisplayLocal( const DistSymmInfo& info, bool beforeFact=true, std::string title="Local DistSymmInfo" )

   Visualize the local reordered nonzero structure before or after symbolic 
   factorization (if Qt5 is not available, just the supernodes are printed).

Print
-----

.. cpp:function:: void Print( const Graph& graph, std::string msg="Graph", std::ostream& os=std::cout )
.. cpp:function:: void Print( const DistGraph& graph, std::string msg="DistGraph", std::ostream& os=std::cout )

   Prints the edges in the graph.

.. cpp:function:: void Print( const SparseMatrix<T>& A, std::string msg="SparseMatrix", std::ostream& os=std::cout )
.. cpp:function:: void Print( const DistSparseMatrix<T>& A, std::string msg="DistSparseMatrix", std::ostream& os=std::cout )

   Prints the nonzero triplets in the sparse matrix.

PrintLocal
----------

.. cpp:function:: void PrintLocal( const DistSymmInfo& info, std::ostream& os=std::cout )

   Print the local supernodal structure.
