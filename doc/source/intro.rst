Introduction
************

Overview
========
Clique is an implementation of the multifrontal algorithm for symmetric and 
Hermitian systems which is designed to be scalable (both in terms of wall time
and memory usage) and as easy to use as possible. As of now, Clique supports 
both :math:`LDL^T` and :math:`LDL^H` factorizations, with or without 
intrafrontal Bunch-Kaufman pivoting, with or without selective inversion, and 
in block or non-block form.

Dependencies
============
Clique is built on top of `Elemental <http://libelemental.org>`_
and `ParMETIS <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview/>`_, 
and both are distributed with Clique. Clique keeps all of Elemental's 
dependencies.

License and copyright
=====================
All files in Clique and Elemental (other than METIS and ParMETIS) are made 
available under the 
`New BSD License <http://www.opensource.org/licenses/bsd-license.php>`_.
The vast majority of files contain the following copyright notice::

    Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
    The University of Texas at Austin, and Stanford University
    All rights reserved.
     
    This file is part of Clique and is under the BSD 2-Clause License, 
    which can be found in the LICENSE file in the root directory, or at 
    http://opensource.org/licenses/BSD-2-Clause
 
Again, METIS and ParMETIS are not distributed under the New BSD License. They
use the following license::

    The ParMETIS/METIS package is copyrighted by the Regents of the
    University of Minnesota. It can be freely used for educational and
    research purposes by non-profit institutions and US government
    agencies only. Other organizations are allowed to use ParMETIS/METIS
    only for evaluation purposes, and any further uses will require prior
    approval. The software may not be sold or redistributed without prior
    approval. One may make copies of the software for their use provided
    that the copies, are not sold or distributed, are used under the same
    terms and conditions.
    
    As unestablished research software, this code is provided on an
    ``as is'' basis without warranty of any kind, either expressed or
    implied. The downloading, or executing any part of this software
    constitutes an implicit agreement to these terms. These terms and
    conditions are subject to change at any time without prior notice.

