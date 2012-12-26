#
#  Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
#  The University of Texas at Austin
# 
#  This file is part of Clique and is under the GNU General Public License,
#  which can be found in the LICENSE file in the root directory, or at 
#  <http://www.gnu.org/licenses/>.
#
cmake_minimum_required(VERSION 2.8.5)
project(Clique)

set(Clique_VERSION_MAJOR 0)
set(Clique_VERSION_MINOR 1)

option(CLIQ_TESTS "Build a collection of test executables" OFF)
option(CLIQ_EXAMPLES "Build a few example drivers" OFF)
option(USE_CUSTOM_ALLTOALLV "Avoid MPI_Alltoallv for performance reasons" ON)
option(BARRIER_IN_ALLTOALLV "Barrier before posting non-blocking recvs" OFF)

if(MANUAL_ELEMENTAL)
  include_directories(${ELEMENTAL_INCLUDE_DIR})
  include_directories(${MPI_CXX_INCLUDE_PATH})
else()
  add_subdirectory(external/elemental)
  include_directories(${PROJECT_BINARY_DIR}/external/elemental/include)
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

option(BUILD_PARMETIS "Build the parallel metis library" ON)
if(BUILD_PARMETIS)
  add_subdirectory(external/parmetis)
  set(HAVE_PARMETIS TRUE)

  set(PARMETIS_ROOT ${PROJECT_SOURCE_DIR}/external/parmetis)
  set(METIS_ROOT ${PROJECT_SOURCE_DIR}/external/parmetis/metis)
  set(GKLIB_INCLUDE_DIR ${PROJECT_BINARY_DIR}/external/parmetis/include)

  include_directories(${GKLIB_INCLUDE_DIR})
  add_subdirectory(src/metis)
  add_subdirectory(src/parmetis)

  include_directories(${METIS_ROOT}/include)
  include_directories(${PARMETIS_ROOT}/include)
elseif(MANUAL_PARMETIS)
  # The following variables need to be specified:
  #   PARMETIS_TLS_PATCH: defined iff there is a GKlib configure file
  #   GKLIB_INCLUDE_DIR: if PARMETIS_TLS_PATCH is defined, this should point
  #                      to the directory containing gklib_tls.h
  #   METIS_ROOT: the full path to the root of the MeTiS source tree
  #   PARMETIS_ROOT: the full path to the root of the ParMeTiS source tree
  set(HAVE_PARMETIS TRUE)
  if(PARMETIS_TLS_PATCH)
    include_directories(${GKLIB_INCLUDE_DIR})
  endif()
  add_subdirectory(src/metis)
  add_subdirectory(src/parmetis)

  include_directories(${METIS_ROOT}/include)
  include_directories(${PARMETIS_ROOT}/include)
endif()

# Create the Clique configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/clique/config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/clique/config.h 
  DESTINATION include/clique)

# Create a file which can be included in Makefile's.
# This is meant to be analogous to PETSc's 'conf/petscvariables' file
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/cliqvariables.cmake
  ${PROJECT_BINARY_DIR}/conf/cliqvariables @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/conf/cliqvariables DESTINATION conf)

# Grab all of the .c, .cpp, .h, and .hpp Clique files
file(GLOB_RECURSE CLIQUE_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/*.cpp")
file(GLOB_RECURSE CLIQUE_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp")
set(CLIQUE_SRC "${CLIQUE_CPP};${CLIQUE_HEADERS}")

# The main library
add_library(clique ${LIBRARY_TYPE} ${CLIQUE_SRC})
if(BUILD_PARMETIS)
  target_link_libraries(clique parmetis-addons metis-addons parmetis)
elseif(MANUAL_PARMETIS)
  target_link_libraries(clique parmetis-addons metis-addons ${PARMETIS_LIBS})
endif()
if(MANUAL_ELEMENTAL)
  target_link_libraries(clique ${ELEM_LIBS})
else()
  target_link_libraries(clique elemental)
endif()
install(TARGETS clique DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${CLIQUE_HEADERS})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}")
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
    list(APPEND TESTS Bisection MultiVectorSolve NestedDissection Refactor
      SimpleVectorSolve SimpleMultiVectorSolve VectorSolve)
  endif()

  # Build the tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests")
  foreach(TEST ${TESTS})
    add_executable(tests-${TEST} ${TEST_DIR}/${TEST}.cpp)
    target_link_libraries(tests-${TEST} clique)
    set_target_properties(tests-${TEST} PROPERTIES
      OUTPUT_NAME ${TEST} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
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
  set(EXAMPLES HelmholtzDirichlet2D HelmholtzDirichlet3D 
               HelmholtzPML2D HelmholtzPML3D)

  # Build the examples
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples")
  foreach(EXAMPLE ${EXAMPLES})
    add_executable(examples-${EXAMPLE} ${EXAMPLE_DIR}/${EXAMPLE}.cpp)
    target_link_libraries(examples-${EXAMPLE} clique)
    set_target_properties(examples-${EXAMPLE} PROPERTIES
      OUTPUT_NAME ${EXAMPLE} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-${EXAMPLE} PROPERTIES
        LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS examples-${EXAMPLE} DESTINATION bin/examples)
  endforeach()
endif()

# If Clique is a subproject, then pass some variables to the parent
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(./cmake/clique_sub.cmake)
endif()
