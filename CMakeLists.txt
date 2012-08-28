#
#  Clique: a scalable implementation of the multifrontal method.
#
#  Copyright (C) 2010-2012 Jack Poulson <jack.poulson@gmail.com>
#  Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
#  The University of Texas at Austin
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
cmake_minimum_required(VERSION 2.8.5)
project(Clique)

set(Clique_VERSION_MAJOR 0)
set(Clique_VERSION_MINOR 1)

option(CLIQ_TESTS "Build a collection of test executables" OFF)
option(CLIQ_EXAMPLES "Build a few example drivers" OFF)
option(USE_CUSTOM_ALLTOALLV "Avoid MPI_Alltoallv for performance reasons" ON)
option(BARRIER_IN_ALLTOALLV "Barrier before posting non-blocking recvs" OFF)

# Add Elemental
add_subdirectory(external/elemental)
include_directories(${PROJECT_BINARY_DIR}/external/elemental/include)
include_directories(${MPI_CXX_INCLUDE_PATH})

option(BUILD_PARMETIS "Build the parallel metis library" ON)
if(BUILD_PARMETIS)
  add_subdirectory(external/parmetis)
  set(HAVE_PARMETIS TRUE)

  include_directories(${PROJECT_BINARY_DIR}/external/parmetis/include)
  add_subdirectory(src/metis)
  add_subdirectory(src/parmetis)

  include_directories(${PROJECT_SOURCE_DIR}/external/parmetis/include)
  include_directories(${PROJECT_SOURCE_DIR}/external/parmetis/metis/include)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Create the Clique configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/clique/config.h 
) 
install(FILES ${PROJECT_BINARY_DIR}/include/clique/config.h 
        DESTINATION include/clique)

# Create a file which can be included in Makefile's.
# This is meant to be analogous to PETSc's 'conf/petscvariables' file
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/cliqvariables.cmake
  ${PROJECT_BINARY_DIR}/conf/cliqvariables @ONLY
)
install(FILES ${PROJECT_BINARY_DIR}/conf/cliqvariables
        DESTINATION conf)

# Grab all of the .c, .cpp, .h, and .hpp Clique files
file(GLOB_RECURSE CLIQUE_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/*.cpp")
file(GLOB_RECURSE CLIQUE_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
     "include/*.h" "include/*.hpp")
set(CLIQUE_SRC "${CLIQUE_CPP};${CLIQUE_HEADERS}")

# The main library
add_library(clique ${LIBRARY_TYPE} ${CLIQUE_SRC})
if(BUILD_PARMETIS)
  target_link_libraries(clique elemental parmetis-addons metis-addons parmetis)
else()
  target_link_libraries(clique elemental)
endif()
set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
install(TARGETS clique DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${CLIQUE_HEADERS})
  add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}"
  )
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_clique_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(clique prepare_clique_headers)

# Make sure the Clique headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(CLIQ_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(TESTS FileSpeed NaturalVectorSolve)
  if(HAVE_PARMETIS)
    list(APPEND TESTS Bisection MultiVectorSolve NestedDissection 
                      SimpleVectorSolve SimpleMultiVectorSolve VectorSolve)
  endif()

  # Build the tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests")
  foreach(TEST ${TESTS})
    add_executable(tests-${TEST} ${TEST_DIR}/${TEST}.cpp)
    target_link_libraries(tests-${TEST} clique)
    set_target_properties(tests-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS tests-${TEST} DESTINATION bin/tests)
  endforeach()
endif()

# Build the example drivers if necessary
if(CLIQ_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)
  set(EXAMPLES HelmholtzDirichlet2D HelmholtzDirichlet3D HelmholtzPML3D)

  # Build the examples
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples")
  foreach(EXAMPLE ${EXAMPLES})
    add_executable(examples-${EXAMPLE} ${EXAMPLE_DIR}/${EXAMPLE}.cpp)
    target_link_libraries(examples-${EXAMPLE} clique)
    set_target_properties(examples-${EXAMPLE} PROPERTIES
                          OUTPUT_NAME ${EXAMPLE}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS examples-${EXAMPLE} DESTINATION bin/examples)
  endforeach()
endif()

################################################################################
# Uncomment if including Clique as a subproject in another build system        #
################################################################################
#set(USE_CUSTOM_ALLTOALLV ${USE_CUSTOM_ALLTOALLV} PARENT_SCOPE)
#set(BARRIER_IN_ALLTOALLV ${BARRIER_IN_ALLTOALLV} PARENT_SCOPE)
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
#set(HAVE_OPENMP ${HAVE_OPENMP} PARENT_SCOPE)
#set(HAVE_F90_INTERFACE ${HAVE_F90_INTERFACE} PARENT_SCOPE)
#set(WITHOUT_PMRRR ${WITHOUT_PMRRR} PARENT_SCOPE)
#set(AVOID_COMPLEX_MPI ${AVOID_COMPLEX_MPI} PARENT_SCOPE)
#set(HAVE_REDUCE_SCATTER_BLOCK ${HAVE_REDUCE_SCATTER_BLOCK} PARENT_SCOPE)
#set(REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE ${REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE} PARENT_SCOPE)
#set(USE_BYTE_ALLGATHERS ${USE_BYTE_ALLGATHERS} PARENT_SCOPE)
