set(METIS_ROOT ${PROJECT_SOURCE_DIR}/external/parmetis/metis)
include_directories(${METIS_ROOT}/GKlib)
include_directories(${METIS_ROOT}/include)
include_directories(${METIS_ROOT}/libmetis)

add_library(metis-addons ${LIBRARY_TYPE} "./Bisect.c")
install(TARGETS metis-addons DESTINATION lib)
