/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CLIQUE_IMPORTS_HPP
#define CLIQUE_IMPORTS_HPP 1

#if defined(BLAS_POST)
#define BLAS(name) name ## _
#else
#define BLAS(name) name
#endif

#if defined(LAPACK_POST)
#define LAPACK(name) name ## _
#else
#define LAPACK(name) name
#endif

#include "clique/imports/blas.hpp"
#include "clique/imports/lapack.hpp"
#include "clique/imports/mpi.hpp"

#endif /* CLIQUE_IMPORTS_HPP */
