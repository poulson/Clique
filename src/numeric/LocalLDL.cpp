/*
   Modification of include/elemental/advanced/LDL.hpp from Elemental.
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

template<typename F> // F represents a real or complex field
void clique::numeric::LocalLDL
( Orientation orientation,
  const symbolic::LocalFactStruct& SLocal,
  const numeric::LocalOrigMatrix<F>& ALocal,
        numeric::LocalFactMatrix<F>& LLocal )
{
#ifndef RELEASE
    PushCallStack("numeric::LocalLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = SLocal.lowerStructs.size();

    // Perform the local factorization
    for( int k=0; k<numSupernodes; ++k )
    {
        
    }
#ifndef RELEASE
    PopCallStack();
#endif
}