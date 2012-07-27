Environment
===========
**TODO**

Set up and clean up
-------------------

.. cpp:function:: void Initialize( int& argc, char**& argv )

   Initialized Clique and (if necessary) Elemental. The usage is very similar
   to ``MPI_Init``, but the command-line information, ``argc`` and ``argv``, 
   do not have to be passed in as pointers:

   .. code-block:: cpp

      #include "clique.hpp"

      int
      main( int argc, char* argv[] )
      {
          cliq::Initialize( argc, argv );

          // ...

          cliq::Finalize();
      }

.. cpp:function:: void Finalize()

   Cleans up all allocated resources from Clique's environment and 
   (if necessary) from Elemental's

.. cpp:function:: bool Initialized()

   Returns whether or not Clique is currently initialized.

Functionality from Elemental
----------------------------

Clique heavily relies on many of 
`Elemental's <http://poulson.github.com/Elemental>`__ data structures, 
particularly the :cpp:class:`DistMatrix\<T,U,V>` class for managing
distributed matrices and vectors. 

.. rubric:: Complex data

.. cpp:type:: struct Complex<R>

   Import of Elemental's `Complex <http://poulson.github.com/Elemental/core/environment.html#Complex:R:>`__ class.

.. cpp:type:: struct Base<F>

   .. cpp:type:: type

      Underlying real datatype of the field ``F``.

.. cpp:function:: typename Base<F>::type Abs( const F& alpha )

   Returns the absolute value of the real or complex variable :math:`\alpha`.

.. cpp:function:: F Sqrt( const F& alpha )

   Returns the square-root of the real or complex variable :math:`\alpha`.

.. rubric:: Classes

.. cpp:class:: Matrix<T>

   Import of Elemental's 
   `Matrix <http://poulson.github.com/Elemental/core/matrix.html>`__ class.

.. cpp:class:: Grid

   Import of Elemental's 
   `Grid <http://poulson.github.com/Elemental/core/grid.html>`__ class.

.. cpp:class:: DistMatrix<T,U,V>

   Import of Elemental's 
   `DistMatrix <http://poulson.github.com/Elemental/core/dist_matrix.html>`__ 
   class.

.. rubric:: Enums

.. cpp:type:: enum Distribution

   See `the Elemental documentation <http://poulson.github.com/Elemental/core/environment.html#Distribution__enum>`__.

.. cpp:type:: enum LeftOrRight

   See `the Elemental documentation <http://poulson.github.com/Elemental/core/environment.html#LeftOrRight__enum>`__.

.. cpp:type:: enum UnitOrNonUnit

   See `the Elemental documentation <http://poulson.github.com/Elemental/core/environment.html#UnitOrNonUnit__enum>`__.

.. cpp:type:: enum Orientation

   See `the Elemental documentation <http://poulson.github.com/Elemental/core/environment.html#Orientation__enum>`__.

.. cpp:type:: enum UpperOrLower

   See `the Elemental documentation <http://poulson.github.com/Elemental/core/environment.html#UpperOrLower__enum>`__.

.. rubric:: Indexing routines

.. cpp:function:: Int Shift( Int rank, Int firstRank, Int numProcs )

   Given a element-wise cyclic distribution over ``numProcs`` processes,
   where the first entry is owned by the process with rank ``firstRank``,
   this routine returns the first entry owned by the process with rank
   ``rank``.

.. cpp:function:: Int LocalLength( Int n, Int shift, Int numProcs )

   Given a vector with :math:`n` entries distributed over ``numProcs``
   processes with shift as defined above, this routine returns the number of
   entries of the vector which are owned by this process.

.. cpp:function:: Int LocalLength( Int n, Int rank, Int firstRank, Int numProcs )

   Given a vector with :math:`n` entries distributed over ``numProcs``
   processes, with the first entry owned by process ``firstRank``, this routine
   returns the number of entries locally owned by the process with rank
   ``rank``.
