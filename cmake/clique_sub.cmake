include(${ELEM_CMAKE_DIR}/elemental_sub.cmake)

set(ELEM_REV ${ELEM_REV} PARENT_SCOPE)
set(ELEM_CMAKE_DIR ${ELEM_CMAKE_DIR} PARENT_SCOPE)
set(USE_CUSTOM_ALLTOALLV ${USE_CUSTOM_ALLTOALLV} PARENT_SCOPE)
set(BARRIER_IN_ALLTOALLV ${BARRIER_IN_ALLTOALLV} PARENT_SCOPE)
set(HAVE_PARMETIS ${HAVE_PARMETIS} PARENT_SCOPE)