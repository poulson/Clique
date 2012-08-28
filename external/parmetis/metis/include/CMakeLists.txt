if(METIS_INSTALL)
  install(FILES metis.h DESTINATION include)
  install(FILES "${PROJECT_BINARY_DIR}/include/gklib_tls.h" DESTINATION include)
endif()
