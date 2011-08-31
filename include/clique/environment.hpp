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
#ifndef CLIQUE_ENVIRONMENT_HPP
#define CLIQUEL_ENVIRONMENT_HPP 1

#include "mpi.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <vector>

#include "clique/config.h"

namespace clique {

typedef unsigned char byte;
 
typedef std::complex<float> scomplex;
typedef std::complex<double> dcomplex;

void Initialize( int& argc, char**& argv );
void Finalize();

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif

}

#include "clique/imports.hpp"

namespace clique {

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename Z> Z Abs( Z alpha );
template<typename Z> Z Abs( std::complex<Z> alpha );
template<typename Z> Z Conj( Z alpha );
template<typename Z> std::complex<Z> Conj( std::complex<Z> alpha );

// We define an output stream that does nothing. This is done so that the 
// root process can be used to print data to a file's ostream while all other 
// processes use a null ostream. 
struct NullStream : std::ostream
{
    struct NullStreamBuffer : std::streambuf
    {
        int overflow( int c ) { return traits_type::not_eof(c); }
    } _nullStreamBuffer;

    NullStream() 
    : std::ios(&_nullStreamBuffer), std::ostream(&_nullStreamBuffer) 
    { }
};

} // clique

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename Z>
inline Z 
clique::Abs( Z alpha )
{ return std::abs(alpha); }

template<typename Z>
inline Z
clique::Abs( std::complex<Z> alpha )
{ return std::abs( alpha ); }

template<typename Z>
inline Z
clique::Conj( Z alpha )
{ return alpha; }

template<typename Z>
inline std::complex<Z>
clique::Conj( std::complex<Z> alpha )
{ return std::conj( alpha ); }

#endif /* CLIQUE_ENVIRONMENT_HPP */

