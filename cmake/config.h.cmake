/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
*/
#ifndef CLIQUE_CONFIG_H
#define CLIQUE_CONFIG_H

/* Version information */
#define Clique_VERSION_MAJOR "@Clique_VERSION_MAJOR@"
#define Clique_VERSION_MINOR "@Clique_VERSION_MINOR@"

/* Miscellaneous */
#cmakedefine USE_CUSTOM_ALLTOALLV
#cmakedefine BARRIER_IN_ALLTOALLV
#cmakedefine HAVE_PARMETIS

#endif /* CLIQUE_CONFIG_H */
