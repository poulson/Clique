/*
   Modification of include/elemental/basic/level3/Trsm/TrsmLLN.hpp 
   from Elemental.
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   Clique: a scalable implementation of the multifrontal algorithm

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
#include "clique.hpp"
using namespace elemental;

template<typename F>
void clique::numeric::SupernodeSolve
( Shape shape, Diagonal diagonal,
  DistMatrix<F,VC,STAR>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::SupernodeSolve");
#endif
    if( shape == LOWER )
        SupernodeLowerSolve( diagonal, A, supernodeSize );
    else
        SupernodeUpperSolve( diagonal, A, supernodeSize );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::SupernodeLowerSolve
( Diagonal diagonal, DistMatrix<F,VC,STAR>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::SupernodeLowerSolve");
#endif
    // HERE
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::SupernodeUpperSolve
( Diagonal diagonal, DistMatrix<F,VC,STAR>& A, int supernodeSize )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::SupernodeUpperSolve");
#endif
    // HERE
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::SupernodeSolve
( Shape shape, Diagonal diagonal,
  DistMatrix<float,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeLowerSolve
( Diagonal diagonal, 
  DistMatrix<float,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeUpperSolve
( Diagonal diagonal, 
  DistMatrix<float,VC,STAR>& A, int supernodeSize );

template void clique::numeric::SupernodeSolve
( Shape shape, Diagonal diagonal,
  DistMatrix<double,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeLowerSolve
( Diagonal diagonal,
  DistMatrix<double,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeUpperSolve
( Diagonal diagonal,
  DistMatrix<double,VC,STAR>& A, int supernodeSize );

template void clique::numeric::SupernodeSolve
( Shape shape, Diagonal diagonal,
  DistMatrix<std::complex<float>,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeLowerSolve
( Diagonal diagonal,
  DistMatrix<std::complex<float>,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeUpperSolve
( Diagonal diagonal,
  DistMatrix<std::complex<float>,VC,STAR>& A, int supernodeSize );

template void clique::numeric::SupernodeSolve
( Shape shape, Diagonal diagonal,
  DistMatrix<std::complex<double>,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeLowerSolve
( Diagonal diagonal,
  DistMatrix<std::complex<double>,VC,STAR>& A, int supernodeSize );
template void clique::numeric::SupernodeUpperSolve
( Diagonal diagonal,
  DistMatrix<std::complex<double>,VC,STAR>& A, int supernodeSize );
