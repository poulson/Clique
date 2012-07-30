The DistMap class
=================
The :cpp:type:`DistMap` class is essentially used for handling distributed 
permutations and can serve as a black box for the vast majority of users.

.. cpp:type:: class DistMap

   .. rubric:: Constructors

   .. cpp:function:: DistMap()

      Constructs a map with zero sources over ``mpi::COMM_WORLD``.

   .. cpp:function:: DistMap( mpi::Comm comm )

      Constructs a map with zero sources over the specified communicator.

   .. cpp:function:: DistMap( int numSources, mpi::Comm comm )

      Constructs a map with the given number of sources over a particular 
      communicator.

   .. rubric:: Map manipulation routines

   .. cpp:function:: void Translate( std::vector<int>& localIndices ) const

      Collectively maps each process's set of local indices.

   .. cpp:function:: void FormInverse( DistMap& inverseMap ) const

      Forms the inverse map.

   .. cpp:function:: void Extend( DistMap& firstMap ) const

      Overwrites the input map with its composition with this map.

   .. cpp:function:: void Extend( const DistMap& firstMap, DistMap& compositeMap ) const

      Sets the composite map equal to `firstMap` composed with this map.

   .. rubric:: High-level information

   .. cpp:function:: int NumSources() const

      The number of sources in the map.

   .. rubric:: Communicator-management

   .. cpp:function:: void SetComm( mpi::Comm comm )

      Reconfigures the map to be distributed over the specified communicator.

   .. cpp:function:: mpi::Comm Comm() const

      The communciator for the distributed map.

   .. rubric:: Distribution information

   .. cpp:function:: int Blocksize() const

      The distribution blocksize of the map: the process with rank ``r``'s
      first local source is global source ``r*blocksize``.

   .. cpp:function:: int FirstLocalSource() const

      The global index of the first source assigned to this process.

   .. cpp:function:: int NumLocalSources() const

      The number of sources assigned to this process.

   .. rubric:: Local data

   .. cpp:function:: T GetLocal( int localSource ) const

      The mapped value of the specified local source.

   .. cpp:function:: void SetLocal( int localSource, int target )
     
      Modifies the mapped value of the specified local source.

   .. cpp:function:: int* LocalBuffer()

      Returns a direct pointer to the mapped values.

   .. cpp:function:: const int* LocalBuffer() const

      Returns a const pointer to the mapped values.

   .. cpp:function:: std::vector<int>& LocalMap()

      Returns the underlying vector of mapped values.

   .. cpp:function:: const std::vector<int>& LocalMap() const

      Returns a const reference to the underlying vector of mapped values.

   .. rubric:: For modifying the size of the map

   .. cpp:function:: void Empty()

      Frees all resources and makes the map have zero sources.

   .. cpp:function:: void ResizeTo( int numSources )

      Resizes the map to have the specified number of sources.

   .. cpp:function:: const DistMap<T>& operator=( const DistMap<T>& x )

      Makes this map a copy of the given map.
