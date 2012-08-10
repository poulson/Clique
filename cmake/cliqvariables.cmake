# To help simplify including Clique in external projects
include @CMAKE_INSTALL_PREFIX@/conf/elemvariables

CLIQ_COMPILE_FLAGS = ${ELEM_COMPILE_FLAGS}
CLIQ_LINK_FLAGS = ${ELEM_LINK_FLAGS}

HAVE_PARMETIS = @HAVE_PARMETIS@
ifeq (${HAVE_PARMETIS},TRUE)
  CLIQ_LIBS = -lclique -lparmetis -lmetis ${ELEM_LIBS}
else
  CLIQ_LIBS = -lclique ${ELEM_LIBS}
endif
