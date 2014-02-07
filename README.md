# Clique


A C++ sparse-direct solver for distributed-memory architectures with an 
emphasis on scalable triangular solves and memory usage.
The library was initially released as part of 
[A parallel sweeping preconditioner for heterogeneous 3D Helmholtz equations](http://epubs.siam.org/doi/abs/10.1137/120871985) ([a preprint is available on arXiv](http://arxiv.org/abs/1204.0111)) and is built directly on top of 
[Elemental](http://github.com/poulson/Elemental).

You may check out a copy of Clique using the command:

    git clone --recursive git://github.com/poulson/Clique.git

The ``--recursive`` argument handles checking out [Elemental](http://github.com/poulson/Elemental) as a [git submodule](http://git-scm.com/book/en/Git-Tools-Submodules). If you have an especially outdated version of git, then you may 
instead need to run:

    git clone git://github.com/poulson/Clique.git
    cd Clique
    git submodule update --init

and, with even older versions of git:

    git clone git://github.com/poulson/Clique.git
    cd Clique
    git submodule init
    git submodule update

### Documentation

The [documentation for the development version of Clique](http://poulson.github.com/Clique) is built using [Sphinx](http://sphinx.pocoo.org).

### Related open-source packages

Implementations:

1. [DSCPACK](http://www.cse.psu.edu/~raghavan/Dscpack/)
2. [MUMPS](http://graal.ens-lyon.fr/MUMPS/)
3. [PSPASES](http://www-users.cs.umn.edu/~mjoshi/pspases/)
4. [SuperLU_Dist](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)

Wrappers:

1. [PETSc](https://www.mcs.anl.gov/petsc/)
2. [Trilinos (Amesos)](http://trilinos.sandia.gov/packages/amesos/)

Note that [PETSc](https://www.mcs.anl.gov/petsc/) contains interfaces for both 
[Clique](http://github.com/poulson/Clique.git) and 
[Elemental](http://github.com/poulson/Elemental.git).
