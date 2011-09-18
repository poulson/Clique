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

struct LocalSymmOrigSupernode
{
    int size, offset;
    int parent; // -1 if root separator
    std::vector<int> children;
    std::vector<int> lowerStruct;
};

struct LocalSymmOrig
{
    std::vector<LocalSymmOrigSupernode> supernodes;
};

struct LocalSymmFactSupernode
{
    bool isLeftChild;
    int size, offset;
    int parent; // -1 if root separator
    std::vector<int> children;

    std::vector<int> lowerStruct;
    std::vector<int> leftChildRelIndices, rightChildRelIndices;
    std::map<int,int> origLowerRelIndices;
};

struct LocalSymmFact
{
    std::vector<LocalSymmFactSupernode> supernodes;
};

struct DistSymmOrigSupernode
{
    int size, offset;
    std::vector<int> lowerStruct;
};

struct DistSymmOrig
{
    mpi::Comm comm;
    std::vector<DistSymmOrigSupernode> supernodes;
};

struct DistSymmFactSupernode
{
    mpi::Comm comm;
    int gridHeight;

    int size, offset;
    std::vector<int> lowerStruct;
    std::map<int,int> origLowerRelIndices;
    std::vector<int> leftChildRelIndices, rightChildRelIndices;

    std::deque<int> leftChildColIndices, leftChildRowIndices,
                    rightChildColIndices, rightChildRowIndices;
    std::vector<int> numChildSendIndices;

    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::deque<int> > childRecvIndices;
};

struct DistSymmFact
{
    std::vector<DistSymmFactSupernode> supernodes;
};

void SymmetricFactorization
( const LocalSymmOrig&   localOrig,
  const DistSymmOrig&    distOrig,
        LocalSymmFact&   localFact,
        DistSymmFact&    distFact, 
        bool storeRecvIndices=true );

void LocalSymmetricFactorization
( const LocalSymmOrig& localOrig,
        LocalSymmFact& localFact );

void DistSymmetricFactorization
( const DistSymmOrig&  distOrig, 
  const LocalSymmFact& localFact, 
        DistSymmFact&  distFact, 
        bool storeRecvIndices=true );

void ComputeRecvIndices( const DistSymmFactSupernode& supernode );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline void SymmetricFactorization
( const LocalSymmOrig& localOrig,
  const DistSymmOrig&  distOrig,
        LocalSymmFact& localFact,
        DistSymmFact&  distFact,
        bool storeRecvIndices )
{
#ifndef RELEASE
    PushCallStack("symbolic::SymmetricFactorization");
#endif
    LocalSymmetricFactorization( localOrig, localFact );
    DistSymmetricFactorization
    ( distOrig, localFact, distFact, storeRecvIndices );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace symbolic
} // namespace clique

#endif /* CLIQUE_SYMBOLIC_SYMM_FACT_HPP */

