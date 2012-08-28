# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckCSourceRuns)

# Setup options.
option(GDB "enable use of GDB" OFF)
option(ASSERT "turn asserts on" OFF)
option(ASSERT2 "additional assertions" OFF)
option(DEBUG "add debugging support" OFF)
option(GPROF "add gprof support" OFF)
option(OPENMP "enable OpenMP support" OFF)
option(PCRE "enable PCRE support" OFF)
option(GKREGEX "enable GKREGEX support" OFF)
option(GKRAND "enable GKRAND support" OFF)

# Add compiler flags.
if(MSVC)
  set(GKlib_COPTS "/Ox")
  set(GKlib_COPTIONS "-DWIN32 -DMSC -D_CRT_SECURE_NO_DEPRECATE -DUSE_GKREGEX")
elseif(MINGW)
  set(GKlib_COPTS "-DUSE_GKREGEX")
else()
  set(GKlib_COPTS "-O3")
  set(GKlib_COPTIONS "-DLINUX -D_FILE_OFFSET_BITS=64")
endif(MSVC)
if(CYGWIN)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DCYGWIN")
endif(CYGWIN)
if(CMAKE_COMPILER_IS_GNUCC)
# GCC opts.
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -std=c99 -fno-strict-aliasing")
  if(NOT MINGW)
      set(GKlib_COPTIONS "${GKlib_COPTIONS} -fPIC")
  endif(NOT MINGW)
# GCC warnings.
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -Wall -pedantic -Wno-unused-variable -Wno-unknown-pragmas")
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Sun")
# Sun insists on -xc99.
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -xc99")
endif(CMAKE_COMPILER_IS_GNUCC)

# Find OpenMP if it is requested.
if(OPENMP)
  include(FindOpenMP)
  if(OPENMP_FOUND)
    set(GKlib_COPTIONS "${GKlib_COPTIONS} -D__OPENMP__ ${OpenMP_C_FLAGS}")
  else()
    message(WARNING "OpenMP was requested but support was not found")
  endif(OPENMP_FOUND)
endif(OPENMP)


# Add various definitions.
if(GDB)
  set(GKlib_COPTS "${GKlib_COPTS} -g")
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -Werror")
endif(GDB)


if(DEBUG)
  set(GKlib_COPTS "-g")
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DDEBUG")
endif(DEBUG)

if(GPROF)
  set(GKlib_COPTS "-pg")
endif(GPROF)

if(NOT ASSERT)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DNDEBUG")
endif(NOT ASSERT)

if(NOT ASSERT2)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DNDEBUG2")
endif(NOT ASSERT2)


# Add various options
if(PCRE)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -D__WITHPCRE__")
endif(PCRE)

if(GKREGEX)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DUSE_GKREGEX")
endif(GKREGEX)

if(GKRAND)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DUSE_GKRAND")
endif(GKRAND)


# Check for features.
check_include_file(execinfo.h HAVE_EXECINFO_H)
if(HAVE_EXECINFO_H)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DHAVE_EXECINFO_H")
endif(HAVE_EXECINFO_H)

check_function_exists(getline HAVE_GETLINE)
if(HAVE_GETLINE)
  set(GKlib_COPTIONS "${GKlib_COPTIONS} -DHAVE_GETLINE")
endif(HAVE_GETLINE)

set(GKlib_TLS "")
# Custom check for TLS
foreach(tls_name "__thread" "__declspec(thread)")
  check_c_source_runs("${tls_name} int x; int main() {x=0; return x;}" HAVE${tls_name})
  if(HAVE${tls_name})
    set(GKlib_TLS "${tls_name}")
    message(STATUS "checking for ${tls_name} thread-local storage - found")
    break()
  else()
    message(STATUS "checking for ${tls_name} thread-local storage - not found")
  endif(HAVE${tls_name})
endforeach()

configure_file(${GKLIB_PATH}/gklib_tls.h.in include/gklib_tls.h)

# Finally set the official C flags.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GKlib_COPTIONS} ${GKlib_COPTS}")

# Find GKlib sources.
file(GLOB GKlib_sources ${GKLIB_PATH}/*.c)
file(GLOB GKlib_includes ${GKLIB_PATH}/*.h)
