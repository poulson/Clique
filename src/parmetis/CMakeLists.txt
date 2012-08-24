set(PARMETIS_ROOT ${PROJECT_SOURCE_DIR}/external/parmetis)
include_directories(${PARMETIS_ROOT}/include)
include_directories(${PARMETIS_ROOT}/libparmetis)
include_directories(${PARMETIS_ROOT}/metis/GKlib)
include_directories(${PARMETIS_ROOT}/metis/include)

add_library(parmetis-addons ${LIBRARY_TYPE} "./ParallelBisect.c")
install(TARGETS parmetis-addons DESTINATION lib)
