/*
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

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDLForwardSolve
( const symbolic::LocalSymmFact& localS,
  const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<F>& localL,
  const numeric::DistSymmFact<F>& distL,
        Matrix<F>& localX )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
#endif
    // HERE
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDLForwardSolve
( const symbolic::LocalSymmFact& localS,
  const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<float>& localL,
  const numeric::DistSymmFact<float>& distL,
        Matrix<float>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::LocalSymmFact& localS,
  const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<double>& localL,
  const numeric::DistSymmFact<double>& distL,
        Matrix<double>& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::LocalSymmFact& localS,
  const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
  const numeric::DistSymmFact<std::complex<float> >& distL,
        Matrix<std::complex<float> >& localX );

template void clique::numeric::DistLDLForwardSolve
( const symbolic::LocalSymmFact& localS,
  const symbolic::DistSymmFact& distS,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
  const numeric::DistSymmFact<std::complex<double> >& distL,
        Matrix<std::complex<double> >& localX );
