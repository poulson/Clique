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
void clique::numeric::LocalLDLForwardSolve
(       symbolic::LocalFactStruct& SLocal, // can't be const due to map...
  const numeric::LocalFactMatrix<F>& LLocal,
  F alpha, Matrix<F>& XLocal )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDLForwardSolve");
#endif
    const int numSupernodes = SLocal.lowerStructs.size();

    for( int k=0; k<numSupernodes; ++k )
    {
        // ???

        // Call the custom supernode forward solve
        // LocalSupernodeLDLForwardSolve
        // ( supernodeSize, alpha, front, XSub );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

