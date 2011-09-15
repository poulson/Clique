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
#ifndef CLIQUE_SYMBOLIC_SYMM_FACT_HPP
#define CLIQUE_SYMBOLIC_SYMM_FACT_HPP 1

namespace clique {
namespace symbolic {

struct LocalOrigStruct
{
    std::vector<int> sizes, offsets;
    std::vector<std::vector<int> > lowerStructs;
    std::vector<std::vector<int> > children;
    std::vector<int> parents; // -1 if root separator
};

struct LocalFactStruct
{
    std::vector<int> sizes, offsets;
    std::vector<std::vector<int> > lowerStructs;
    std::vector<std::vector<int> > children;
    std::vector<int> parents; // -1 if root separator

    std::vector<std::map<int,int> > origLowerRelIndices;
    std::vector<std::vector<int> > leftChildRelIndices, rightChildRelIndices;
};

struct DistOrigStruct
{
    mpi::Comm comm;
    std::vector<int> sizes, offsets;
    std::vector<std::vector<int> > lowerStructs;
};

struct DistFactStruct
{
    std::vector<mpi::Comm> comms;

    std::vector<int> sizes, offsets;
    std::vector<std::vector<int> > lowerStructs;
    std::vector<std::map<int,int> > origLowerRelIndices;
    std::vector<std::vector<int> > leftChildRelIndices, rightChildRelIndices;
};

void SymmetricFactorization
( const LocalOrigStruct& localOrig,
  const DistOrigStruct&  distOrig,
        LocalFactStruct& localFact,
        DistFactStruct&  distFact );

void LocalSymmetricFactorization
( const LocalOrigStruct& localOrig,
        LocalFactStruct& localFact );

void DistSymmetricFactorization
( const DistOrigStruct&  distOrig, 
  const LocalFactStruct& localFact, 
        DistFactStruct&  distFact );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline void SymmetricFactorization
( const LocalOrigStruct& localOrig,
  const DistOrigStruct&  distOrig,
        LocalFactStruct& localFact,
        DistFactStruct&  distFact )
{
#ifndef RELEASE
    PushCallStack("symbolic::SymmetricFactorization");
#endif
    LocalSymmetricFactorization( localOrig, localFact );
    DistSymmetricFactorization( distOrig, localFact, distFact );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace symbolic
} // namespace clique

#endif /* CLIQUE_SYMBOLIC_SYMM_FACT_HPP */

