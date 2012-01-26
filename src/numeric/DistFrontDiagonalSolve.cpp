/*
   Modification of include/elemental/basic/level3/Trsm/TrsmLLN.hpp 
   from Elemental.
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
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

namespace cliq {

template<typename F>
void numeric::DistFrontDiagonalSolve
( const DistMatrix<F,VC,STAR>& d, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontDiagonalSolve");
    if( d.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("Incompatible alignments");
    if( d.Width() != 1 )
        throw std::logic_error("d must be a column vector");
    if( d.Height() != X.Height() )
        throw std::logic_error("Invalid height of X");
#endif
    elem::DiagonalSolve
    ( LEFT, NORMAL, d.LockedLocalMatrix(), X.LocalMatrix(), true );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

template void cliq::numeric::DistFrontDiagonalSolve
( const DistMatrix<float,VC,STAR>& d, 
        DistMatrix<float,VC,STAR>& X );

template void cliq::numeric::DistFrontDiagonalSolve
( const DistMatrix<double,VC,STAR>& d, 
        DistMatrix<double,VC,STAR>& X );

template void cliq::numeric::DistFrontDiagonalSolve
( const DistMatrix<Complex<float>,VC,STAR>& d, 
        DistMatrix<Complex<float>,VC,STAR>& X );

template void cliq::numeric::DistFrontDiagonalSolve
( const DistMatrix<Complex<double>,VC,STAR>& d, 
        DistMatrix<Complex<double>,VC,STAR>& X );
