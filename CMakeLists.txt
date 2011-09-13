#
#  Clique: a scalable implementation of the multifrontal method.
#
#  Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
#  Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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

option(BUILD_TESTS "Build a collection of test executables" OFF)

# We don't need a Hermitian eigensolver
option(BUILD_PMRRR "Build the parallel tridiagonal eigensolver" OFF)
add_subdirectory(external/elemental)
include_directories("${PROJECT_BINARY_DIR}/external/elemental/include")
include_directories(${MPI_CXX_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Create the Clique configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/include/clique/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/clique/config.h 
) 
install(FILES ${PROJECT_BINARY_DIR}/include/clique/config.h 
        DESTINATION include/clique)

# Grab all of the .c, .cpp, .h, and .hpp Clique files
file(GLOB_RECURSE CLIQUE_CPP RELATIVE ${PROJECT_SOURCE_DIR} 
     "src/*.c" "src/*.cpp")
file(GLOB_RECURSE CLIQUE_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
     "include/*.h" "include/*.hpp")
set(CLIQUE_SRC "${CLIQUE_CPP};${CLIQUE_HEADERS}")

# The main library
add_library(clique STATIC ${CLIQUE_SRC})
target_link_libraries(clique elemental)
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
endforeach(HEADER)
add_custom_target(prepare_clique_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(clique prepare_clique_headers)

# Make sure the Clique headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(CORE_TESTS )

  # Build the core tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/core")
  foreach(TEST ${CORE_TESTS})
    add_executable(tests-core-${TEST} ${TEST_DIR}/${TEST}.cpp)
    target_link_libraries(tests-core-${TEST} clique cmake-dummy-lib)
    set_target_properties(tests-core-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-core-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-core-${TEST} DESTINATION bin/tests/core)
  endforeach(TEST)
endif(BUILD_TESTS)
