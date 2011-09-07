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

if(CMAKE_BUILD_TYPE MATCHES Release)
  option(RELEASE "Build high-performance version of library" ON)
endif(CMAKE_BUILD_TYPE MATCHES Release)

option(BUILD_TESTS "Build a collection of test executables" OFF)
option(AVOID_COMPLEX_MPI "Avoid complex MPI routines for robustness" ON)
option(USE_BYTE_ALLGATHERS "Avoid a known MPI-performance bug on BG/P" ON)

find_package(MPI)
if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif(NOT MPI_CXX_FOUND)
include_directories(${MPI_CXX_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILER_FLAGS}")

# Query the size of a void pointer in order to determine whether or not this is
# a 32-bit system.
# Unfortunately, CMAKE_SIZEOF_VOID_P seems to be incorrect on some Windows 
# platforms, so I am switching to directly checking the size of void*
include(CheckTypeSize)
if(NOT SIZEOF_VOID_P)
  check_type_size("void*" SIZEOF_VOID_P)
  if(${SIZEOF_VOID_P} EQUAL 4)
    set(32_BIT_SYSTEM TRUE)
  else(${SIZEOF_VOID_P} EQUAL 4)
    set(32_BIT_SYSTEM FALSE)    
  endif(${SIZEOF_VOID_P} EQUAL 4)
endif(NOT SIZEOF_VOID_P)

# Since the CMake find_blas for some reason requires a Fortran compiler, we
# must build our own.
if(MATH_LIBS)
  set(NEED_MATH FALSE)
  message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
  set(NEED_MATH TRUE)
endif(MATH_LIBS)
if(NEED_MATH)
  if(PURE)
    set(MATH_DESC "Unthreaded BLAS/LAPACK link flags")
  else(PURE)
    set(MATH_DESC "Threaded BLAS/LAPACK link flags")
  endif(PURE)
  if(APPLE)
    if(NEED_MATH)
      set(MATH_LIBS "-framework vecLib" CACHE STRING ${MATH_DESC})
      message(STATUS "Using Apple vecLib framework.")
    endif(NEED_MATH)
  else(APPLE)
    set(FOUND_MATH_LIBS FALSE)
    # Look for MKL first
    if(32_BIT_SYSTEM)
      set(INTEL_ARCH_SUBDIRS 32/lib lib/32)
      set(MKL_ARCH_SUBDIRS 32/lib lib/32)
    else(32_BIT_SYSTEM)
      set(INTEL_ARCH_SUBDIRS em64t/lib intel64/lib lib/em64t lib/intel64)
      set(MKL_ARCH_SUBDIRS em64t/lib intel64/lib lib/em64t lib/intel64 
          mkl/em64t/lib mkl/intel64/lib mkl/lib/em64t mkl/lib/intel64)
    endif(32_BIT_SYSTEM)
    if(INTEL_ROOT)
      message(STATUS "Will also search INTEL_ROOT=${INTEL_ROOT} for MKL")
      message(STATUS "INTEL_ARCH_SUBDIRS=${INTEL_ARCH_SUBDIRS}")
    endif(INTEL_ROOT)
    if(MKL_ROOT)
      message(STATUS "Will also search MKL_ROOT=${MKL_ROOT} for MKL")
      message(STATUS "MKL_ARCH_SUBDIRS=${INTEL_ARCH_SUBDIRS}")
    endif(MKL_ROOT)
    # There is a cyclic dependency between MKL_CORE and 
    # MKL_SEQUENTIAL/MKL_INTEL_THREAD with some versions of MKL
    if(NOT WIN32)
      find_library(PTHREAD_LIB pthread DOC "Pthreads library")
    endif(NOT WIN32)
    find_library(
      CMATH_LIB
      NAMES m libmmd.lib
      PATHS ${INTEL_ROOT} 
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
      DOC "C math library"
    )
    find_library(
      IOMP5_LIB
      NAMES libiomp5.a libiomp5md.lib
      PATHS ${INTEL_ROOT} 
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    )
    find_library(
      GUIDE_LIB
      NAMES libguide.a libguide.lib
      PATHS ${INTEL_ROOT}
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    )
    find_library(
      MKL_CORE_LIB
      NAMES libmkl_core.a mkl_core.lib 
      PATHS ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    find_library(
      MKL_INTEL_LIB
      NAMES libmkl_intel_lp64.a mkl_intel_lp64.lib 
      PATHS ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    if(PURE)
      find_library(
        MKL_SEQUENTIAL_LIB
        NAMES libmkl_sequential.a mkl_sequential.lib
        PATHS ${MKL_ROOT}
        PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
      )
      set(MKL_REQUIRED MKL_INTEL MKL_SEQUENTIAL MKL_CORE MKL_SEQUENTIAL CMATH
                       GUIDE IOMP5)
    else(PURE)
      find_library(
        MKL_INTEL_THREAD_LIB
        NAMES libmkl_intel_thread.a mkl_intel_thread.lib
        PATHS ${MKL_ROOT}
        PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
      )
      set(MKL_REQUIRED MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_INTEL_THREAD 
                       CMATH GUIDE IOMP5)
    endif(PURE)
    if(NOT WIN32)
      list(APPEND MKL_REQUIRED PTHREAD)
    endif(NOT WIN32)
    set(MKL_FOUND TRUE)
    set(MATH_LIBS "")
    foreach(NAME ${MKL_REQUIRED})
      if(${NAME}_LIB)
        message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
        list(APPEND MATH_LIBS ${${NAME}_LIB})
      else(${NAME}_LIB)
        message(STATUS "Could not find ${NAME}_LIB") 
        set(MATH_LIBS "")
        set(MKL_FOUND FALSE)
      endif(${NAME}_LIB)
    endforeach(NAME)
    if(MKL_FOUND)
      message(STATUS "Using MKL's math libraries.")
      message(STATUS "MATH_LIBS=${MATH_LIBS}")
      set(NEED_MATH FALSE)
    endif(MKL_FOUND)
    if(NEED_MATH)
      # Look for default BLAS and LAPACK
      if(REFERENCE_ROOT)
        message(STATUS 
          "Will search REFERENCE_ROOT=${REFERENCE_ROOT} for math libs")
      endif(REFERENCE_ROOT)
      set(REFERENCE_REQUIRED LAPACK BLAS)
      find_library(BLAS_LIB NAMES blas PATHS ${REFERENCE_ROOT})
      find_library(LAPACK_LIB NAMES lapack reflapack PATHS ${REFERENCE_ROOT}) 
      set(REFERENCE_FOUND TRUE)
      foreach(NAME ${REFERENCE_REQUIRED})
        if(${NAME}_LIB)
          message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
          list(APPEND MATH_LIBS ${${NAME}_LIB})
        else(${NAME}_LIB)
          message(STATUS "Could not find ${NAME}_LIB")
          set(MATH_LIBS "")
          set(REFERENCE_FOUND FALSE)
        endif(${NAME}_LIB)
      endforeach(NAME)
      if(REFERENCE_FOUND)
        message(WARNING "Using reference BLAS/LAPACK.")
      else(REFERENCE_FOUND) 
        message(FATAL_ERROR 
          "Could not find BLAS/LAPACK libraries. Please specify MATH_LIBS")
      endif(REFERENCE_FOUND)
    endif(NEED_MATH)
  endif(APPLE)
endif(NEED_MATH)

# Attempt to detect the BLAS and LAPACK underscore conventions. 
# We currently only handle whether or not there is an underscore appended.
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(daxpy HAVE_DAXPY)
if(HAVE_DAXPY)
  set(BLAS_POST FALSE)
  set(BLAS_DEFS "")
else(HAVE_DAXPY)
  check_function_exists(daxpy_ HAVE_DAXPY_POST)
  if(HAVE_DAXPY_POST)
    set(BLAS_POST TRUE)
    set(BLAS_DEFS "-DBLAS_POST")
  else(HAVE_DAXPY_POST)
    message(FATAL_ERROR "Could not determine BLAS format.")
  endif(HAVE_DAXPY_POST)
endif(HAVE_DAXPY)
check_function_exists(dpotrf HAVE_DPOTRF)
if(HAVE_DPOTRF)
  set(LAPACK_POST FALSE)
  set(LAPACK_DEFS "")
else(HAVE_DPOTRF)
  check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
  if(HAVE_DPOTRF_POST)
    set(LAPACK_POST TRUE)
    set(LAPACK_DEFS "-DLAPACK_POST")
  else(HAVE_DPOTRF_POST)
    message(FATAL_ERROR "Could not determine LAPACK format.")
  endif(HAVE_DPOTRF_POST)
endif(HAVE_DPOTRF)

# Look for restrict support
include(CheckCXXSourceCompiles)
set(RESTRICT_CODE
    "int main(void)
     {
         int* RESTRICT a;
         return 0;
     }
    ")
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
if(HAVE___restrict__)
  set(RESTRICT "__restrict__")
  message(STATUS "Using __restrict__ keyword.")
else(HAVE___restrict__)
  set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict")
  check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict)
  if(HAVE___restrict)
    set(RESTRICT "__restrict")
    message(STATUS "Using __restrict keyword.")
  else(HAVE___restrict)
    set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=restrict")
    check_cxx_source_compiles("${RESTRICT_CODE}" HAVE_restrict)
    if(HAVE_restrict)
      set(RESTRICT "restrict")
      message(STATUS "Using restrict keyword.")
    else(HAVE_restrict)
      set(RESTRICT "")
      message(STATUS "Could not find a restrict keyword.")
    endif(HAVE_restrict)
  endif(HAVE___restrict)
endif(HAVE___restrict__)

# Create the Clique configuration header
configure_file( 
  ${CMAKE_SOURCE_DIR}/include/clique/config.h.cmake
  ${CMAKE_BINARY_DIR}/include/clique/config.h 
) 
install(FILES ${CMAKE_BINARY_DIR}/include/clique/config.h 
        DESTINATION include/clique)

# Grab all of the .c, .cpp, .h, and .hpp Clique files
file(GLOB_RECURSE CLIQUE_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "src/*.c" "src/*.cpp")
file(GLOB_RECURSE CLIQUE_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "include/*.h" "include/*.hpp")
set(CLIQUE_SRC "${CLIQUE_CPP};${CLIQUE_HEADERS}")

# Copy the headers into the build directory
foreach(HEADER ${CLIQUE_HEADERS})
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}
            ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
  )
  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach(HEADER)

# Make sure the Clique headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Create a dummy library in order to be able to force the math (and potentially
# MPI libraries) to be linked last
add_library(cmake-dummy-lib STATIC cmake/CMakeDummyFunction.cpp)
target_link_libraries(cmake-dummy-lib ${MATH_LIBS} ${MPI_CXX_LIBRARIES})

# Build the main library
add_library(clique STATIC ${CLIQUE_SRC})
set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
install(TARGETS clique DESTINATION lib)

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
