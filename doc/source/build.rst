Build system
************
Clique's build system sits on top of 
`Elemental's <http://poulson.github.com/Elemental/build.html>`_, and so it is 
a good idea to familiarize yourself with Elemental's build system and 
dependencies first.

Getting Clique's source
=======================
The best way to get a copy of Clique's source is to install 
`Mercurial <http://mercurial.selenic.com>`_ and run ::

    hg clone http://bitbucket.org/poulson/clique clique

Building Clique
===============
Clique's build works essentially the same as Elemental's, which is described 
`here <http://poulson.github.com/Elemental/build.html#building-elemental>`_.
Assuming all `dependencies <http://poulson.github.com/Elemental/build.html#dependencies>`_ 
have been installed, Clique can often be built and installed using the commands::

    cd clique
    mkdir build
    cd build
    cmake ..
    make
    make install

Note that the default install location is system-wide, e.g., ``/usr/local``.
The installation directory can be changed at any time by running::

    cmake -D CMAKE_INSTALL_PREFIX=/your/desired/install/path ..
    make install

Testing the installation
========================
Once Clique has been installed, it is easy to test whether or not it is 
functioning. Assuming that Clique's source code sits in the directory ``/home/username/clique``, and that Clique was installed in ``/usr/local``, then one can
create the following Makefile from any directory::

    include /usr/local/conf/cliqvariables

    DistSparseMatrix: /home/username/clique/tests/DistSparseMatrix.cpp
        ${CXX} ${CLIQ_COMPILE_FLAGS} $< -o $@ ${CLIQ_LINK_FLAGS} ${CLIQ_LIBS}

and then simply running ``make`` should build the test driver.

You can also build a handful of test drivers by using the CMake option::

    -D CLIQ_TESTS=ON

Clique as a subproject
======================
Adding Clique as a dependency into a project which uses CMake for its build 
system is relatively straightforward: simply put an entire copy of the 
Clique source tree in a subdirectory of your main project folder, say 
``external/clique``, and uncomment out the bottom section of Clique's 
``CMakeLists.txt``, i.e., change ::

    ################################################################################
    # Uncomment if including Clique as a subproject in another build system        #
    ################################################################################
    #set(USE_CUSTOM_ALLTOALLV ${USE_CUSTOM_ALLTOALLV} PARENT_SCOPE)
    #set(HAVE_PARMETIS ${HAVE_PARMETIS} PARENT_SCOPE)
    #
    #set(LIBRARY_TYPE ${LIBRARY_TYPE} PARENT_SCOPE)
    #set(MPI_C_COMPILER ${MPI_C_COMPILER} PARENT_SCOPE)
    #set(MPI_C_INCLUDE_PATH ${MPI_C_INCLUDE_PATH} PARENT_SCOPE)
    #set(MPI_C_COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS} PARENT_SCOPE)
    #set(MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS} PARENT_SCOPE)
    #set(MPI_C_LIBRARIES ${MPI_C_LIBRARIES} PARENT_SCOPE)
    #set(MPI_CXX_COMPILER ${MPI_CXX_COMPILER} PARENT_SCOPE)
    #set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} PARENT_SCOPE)
    #set(MPI_CXX_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS} PARENT_SCOPE)
    #set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} PARENT_SCOPE)
    #set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} PARENT_SCOPE)
    #set(MATH_LIBS ${MATH_LIBS} PARENT_SCOPE)
    #set(RESTRICT ${RESTRICT} PARENT_SCOPE)
    #set(RELEASE ${RELEASE} PARENT_SCOPE)
    #set(BLAS_POST ${BLAS_POST} PARENT_SCOPE)
    #set(LAPACK_POST ${LAPACK_POST} PARENT_SCOPE)
    #set(HAVE_F90_INTERFACE ${HAVE_F90_INTERFACE} PARENT_SCOPE)
    #set(WITHOUT_PMRRR ${WITHOUT_PMRRR} PARENT_SCOPE)
    #set(AVOID_COMPLEX_MPI ${AVOID_COMPLEX_MPI} PARENT_SCOPE)
    #set(HAVE_REDUCE_SCATTER_BLOCK ${HAVE_REDUCE_SCATTER_BLOCK} PARENT_SCOPE)
    #set(REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE ${REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE} PARENT_SCOPE)
    #set(USE_BYTE_ALLGATHERS ${USE_BYTE_ALLGATHERS} PARENT_SCOPE)
    
to ::

    ################################################################################
    # Uncomment if including Clique as a subproject in another build system        #
    ################################################################################
    set(USE_CUSTOM_ALLTOALLV ${USE_CUSTOM_ALLTOALLV} PARENT_SCOPE)
    set(HAVE_PARMETIS ${HAVE_PARMETIS} PARENT_SCOPE)
    
    set(LIBRARY_TYPE ${LIBRARY_TYPE} PARENT_SCOPE)
    set(MPI_C_COMPILER ${MPI_C_COMPILER} PARENT_SCOPE)
    set(MPI_C_INCLUDE_PATH ${MPI_C_INCLUDE_PATH} PARENT_SCOPE)
    set(MPI_C_COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS} PARENT_SCOPE)
    set(MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS} PARENT_SCOPE)
    set(MPI_C_LIBRARIES ${MPI_C_LIBRARIES} PARENT_SCOPE)
    set(MPI_CXX_COMPILER ${MPI_CXX_COMPILER} PARENT_SCOPE)
    set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} PARENT_SCOPE)
    set(MPI_CXX_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS} PARENT_SCOPE)
    set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} PARENT_SCOPE)
    set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} PARENT_SCOPE)
    set(MATH_LIBS ${MATH_LIBS} PARENT_SCOPE)
    set(RESTRICT ${RESTRICT} PARENT_SCOPE)
    set(RELEASE ${RELEASE} PARENT_SCOPE)
    set(BLAS_POST ${BLAS_POST} PARENT_SCOPE)
    set(LAPACK_POST ${LAPACK_POST} PARENT_SCOPE)
    set(HAVE_F90_INTERFACE ${HAVE_F90_INTERFACE} PARENT_SCOPE)
    set(WITHOUT_PMRRR ${WITHOUT_PMRRR} PARENT_SCOPE)
    set(AVOID_COMPLEX_MPI ${AVOID_COMPLEX_MPI} PARENT_SCOPE)
    set(HAVE_REDUCE_SCATTER_BLOCK ${HAVE_REDUCE_SCATTER_BLOCK} PARENT_SCOPE)
    set(REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE ${REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE} PARENT_SCOPE)
    set(USE_BYTE_ALLGATHERS ${USE_BYTE_ALLGATHERS} PARENT_SCOPE)

Afterwards, create a ``CMakeLists.txt`` in your main project folder that builds
off of the following snippet::

    cmake_minimum_required(VERSION 2.8.5)
    project(Foo)

    add_subdirectory(external/clique)
    include_directories("${PROJECT_BINARY_DIR}/external/clique/include")
    if(HAVE_PARMETIS)
      include_directories(
        "${PROJECT_SOURCE_DIR}/external/clique/external/parmetis/include"
      )
      include_directories(
        "${PROJECT_SOURCE_DIR}/external/clique/external/parmetis/metis/include"
      )
    endif()
    include_directories(
      "${PROJECT_BINARY_DIR}/external/clique/external/elemental/include")
    )
    include_directories(${MPI_CXX_INCLUDE_PATH})
     
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    
    # Build your project here
    # e.g.,
    #   add_library(foo STATIC ${FOO_SRC})
    #   target_link_libraries(foo clique)

Troubleshooting
===============
If you run into problems, please email
`jack.poulson@gmail.com <mailto:jack.poulson@gmail.com>`_. If you are having 
build problems, please make sure to attach the file ``include/clique/config.h``,
which should be generated within your build directory.
