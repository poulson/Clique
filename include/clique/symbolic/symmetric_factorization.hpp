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

struct LocalSymmStructure
{
    mpi::Comm comm;
    std::vector<int> sizes, offsets;
    std::vector<std::vector<int> > lowerStructs;
};

struct LocalSymmFact
{
    mpi::Comm comm;
    std::vector<int> sizes, offsets;
    std::vector<std::vector<int> > lowerStructs;
    std::vector<std::vector<int> > leftChildMaps, rightChildMaps;
};

void DistSymmetricFactorization
( const LocalSymmStructure& symmStruct, LocalSymmFact& symmFact );

} // namespace symbolic
} // namespace clique

#endif /* CLIQUE_SYMBOLIC_SYMM_FACT_HPP */
