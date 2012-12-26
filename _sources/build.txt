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

    VectorSolve: /home/username/clique/tests/VectorSolve.cpp
        ${CXX} ${CLIQ_COMPILE_FLAGS} $< -o $@ ${CLIQ_LINK_FLAGS} ${CLIQ_LIBS}

and then simply running ``make`` should build the test driver.

You can also build a handful of test drivers by using the CMake option::

    -D CLIQ_TESTS=ON

Clique as a subproject
======================
Adding Clique as a dependency into a project which uses CMake for its build 
system is relatively straightforward: simply put an entire copy of the 
Clique source tree in a subdirectory of your main project folder, say 
``external/clique``, 
and then create a ``CMakeLists.txt`` in your main project folder that builds
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
     
    # Build your project here
    # e.g.,
    #   add_library(foo ${LIBRARY_TYPE} ${FOO_SRC})
    #   target_link_libraries(foo clique)

Troubleshooting
===============
If you run into problems, please email
`jack.poulson@gmail.com <mailto:jack.poulson@gmail.com>`_. If you are having 
build problems, please make sure to attach the file ``include/clique/config.h``,
which should be generated within your build directory.
