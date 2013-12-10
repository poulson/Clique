Environment
===========

This section describes the routines and data structures which help set up 
Clique's programming environment: it discusses initialization of Clique,
call stack manipulation, and data structures and routines imported from 
`Elemental <https://libelemental.org>`__.

Set up and clean up
-------------------

.. cpp:function:: void Initialize( int& argc, char**& argv )

   Initializes Clique and (if necessary) Elemental. The usage is very similar
   to ``MPI_Init``, but the command-line information, ``argc`` and ``argv``, 
   does not have to be passed in as pointers:

   .. code-block:: cpp

      #include "clique.hpp"

      int
      main( int argc, char* argv[] )
      {
          cliq::Initialize( argc, argv );
          // ...
          cliq::Finalize();
          return 0;
      }

.. cpp:function:: void Finalize()

   Cleans up all allocated resources from Clique's environment and 
   (if necessary) from Elemental's

.. cpp:function:: bool Initialized()

   Returns whether or not Clique is currently initialized.

Call stack manipulation
-----------------------

.. note:: 

   The following call stack manipulation routines are only available in
   non-release builds (i.e., PureDebug and HybridDebug) and are meant to allow
   for the call stack to be printed (via :cpp:func:`DumpCallStack`) when an 
   exception is caught.

.. cpp:type:: class CallStackEntry

   Typically used in the following form:

   .. code-block:: cpp

      void SomeRoutine()
      {
          DEBUG_ONLY(CallStackEntry cse("SomeRoutine"))
          // Routine implemented here
      }

.. cpp:function:: void DumpCallStack()

   Print (and empty) the contents of the call stack.

Functionality from Elemental
----------------------------

Clique heavily relies on many of 
`Elemental's <http://libelemental.org>`__ data structures, 
particularly the :cpp:type:`DistMatrix\<T,U,V>` class for managing
distributed matrices and vectors. 

.. rubric:: Complex data

.. cpp:type:: struct Complex<R>

   This is currently a typedef to C++11's ``std::complex``.

.. cpp:type:: template Base<F>

   Template alias for underlying real datatype of the field ``F``.

.. cpp:function:: Base<F> Abs( const F& alpha )

   Returns the absolute value of the real or complex variable :math:`\alpha`.

.. cpp:function:: F Sqrt( const F& alpha )

   Returns the square-root of the real or complex variable :math:`\alpha`.

.. rubric:: Classes

.. cpp:type:: class Matrix<T>

   Import of Elemental's 
   `Matrix <http://libelemental.org/documentation/dev/core/matrix.html>`__ 
   class.

.. cpp:type:: class Matrix<F>

   Same as above, but this implies that the underlying datatype `F` is a field.

.. cpp:type:: class Grid

   Import of Elemental's 
   `Grid <http://libelemental.org/documentation/dev/core/grid.html>`__ class.

.. cpp:type:: class DistMatrix<T,U,V>

   Import of Elemental's 
   `DistMatrix <http://libelemental.org/documentation/dev/core/dist_matrix.html>`__ 
   class.

.. cpp:type:: class DistMatrix<F,U,V>

   Same as above, but this implies that the underlying datatype `F` is a field.

.. cpp:type:: class DistMatrix<T>

.. cpp:type:: class DistMatrix<T,MC,MR>

   A partial specialization of the :cpp:type:`DistMatrix\<T,U,V>` class to the 
   `standard matrix distribution <http://libelemental.org/documentation/dev/core/dist_matrix.html#mc-mr>`__, ``[MC,MR]``.

.. cpp:type:: class DistMatrix<F>

.. cpp:type:: class DistMatrix<F,MC,MR>

   Same as above, but this implies that the underlying datatype `F` is a field.

.. cpp:type:: class DistMatrix<T,VC,STAR>

   A partial specialization of the :cpp:type:`DistMatrix\<T,U,V>` class to a 
   `column-major vector distribution <http://libelemental.org/documentation/dev/core/dist_matrix.html#vc>`__, ``[VC,*]``.

.. cpp:type:: class DistMatrix<F,VC,STAR>

   Same as above, but this implies that the underlying datatype `F` is a field.

.. rubric:: Imported libraries

Elemental provides 
`high-level interfaces to several libraries <http://libelemental.org/documentation/dev/core/imports.html>`__, 
and several of those interfaces are used within Clique. 

* `BLAS <http://libelemental.org/documentation/dev/core/imports/blas.html>`__: 
  exposed in the ``cliq::blas`` namespace

* `LAPACK <http://libelemental.org/documentation/dev/core/imports/lapack.html>`__:
  exposed in the ``cliq::lapack`` namespace

* `MPI <http://libelemental.org/documentation/dev/core/imports/mpi.html>`__:
  exposed in the ``cliq::mpi`` namespace

.. rubric:: Enums

.. cpp:type:: enum Distribution

   See `the Elemental documentation <http://libelemental.org/documentation/core/environment.html#Distribution__enum>`__.

.. cpp:type:: enum LeftOrRight

   See `the Elemental documentation <http://libelemental.org/documentation/dev/core/environment.html#LeftOrRight__enum>`__.

.. cpp:type:: enum UnitOrNonUnit

   See `the Elemental documentation <http://libelemental.org/documentation/dev/core/environment.html#UnitOrNonUnit__enum>`__.

.. cpp:type:: enum Orientation

   See `the Elemental documentation <http://libelemental.org/documentation/dev/core/environment.html#Orientation__enum>`__.

.. cpp:type:: enum UpperOrLower

   See `the Elemental documentation <http://libelemental.org/documentation/dev/core/environment.html#UpperOrLower__enum>`__.

.. rubric:: Indexing routines

.. cpp:function:: int Shift( int rank, int firstRank, int numProcs )

   Given a element-wise cyclic distribution over ``numProcs`` processes,
   where the first entry is owned by the process with rank ``firstRank``,
   this routine returns the first entry owned by the process with rank
   ``rank``.

.. cpp:function:: int Length( int n, int shift, int numProcs )

   Given a vector with :math:`n` entries distributed over ``numProcs``
   processes with shift as defined above, this routine returns the number of
   entries of the vector which are owned by this process.

.. cpp:function:: int Length( int n, int rank, int firstRank, int numProcs )

   Given a vector with :math:`n` entries distributed over ``numProcs``
   processes, with the first entry owned by process ``firstRank``, this routine
   returns the number of entries locally owned by the process with rank
   ``rank``.
