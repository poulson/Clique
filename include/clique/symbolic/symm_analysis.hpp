/*
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
#ifndef CLIQUE_SYMM_ANALYSIS_HPP
#define CLIQUE_SYMM_ANALYSIS_HPP 1

namespace cliq {

struct LocalSymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, offset; 
    int parent; // -1 if root separator
    std::vector<int> children;
    std::vector<int> origLowerStruct;

    //
    // The following is computed during analysis
    //
    bool isLeftChild;
    int myOffset;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelIndices;
    // The relative indices of the left and right children
    // (maps from the child update indices to our frontal indices).
    std::vector<int> leftChildRelIndices, rightChildRelIndices;
};

struct LocalSymmInfo
{
    std::vector<LocalSymmNodeInfo> nodes;
};

struct DistSymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, offset;
    std::vector<int> origLowerStruct;
    mpi::Comm comm;

    //
    // The following is computed during analysis
    //
    Grid* grid;
    int myOffset, leftChildSize, rightChildSize;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelIndices;

    // The relative indices of the left and right children
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftChildRelIndices, rightChildRelIndices;

    // Helpers for the factorization
    std::vector<int> numChildFactSendIndices;
    std::deque<int> leftChildFactColIndices, leftChildFactRowIndices,
                    rightChildFactColIndices, rightChildFactRowIndices;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::deque<int> > childFactRecvIndices;

    // Helpers for solving with 1d right-hand sides
    std::deque<int> leftChildSolveIndices, rightChildSolveIndices;
    int localSize1d, localOffset1d;
    std::vector<int> numChildSolveSendIndices;
    std::vector<std::deque<int> > childSolveRecvIndices;
};

struct DistSymmInfo
{
    std::vector<DistSymmNodeInfo> nodes;
};

struct SymmInfo
{
    LocalSymmInfo local;
    DistSymmInfo dist;
};

void SymmetricAnalysis
( const SymmElimTree& eTree, SymmInfo& info, 
  bool storeFactRecvIndices=true );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

void LocalSymmetricAnalysis
( const SymmElimTree& eTree, SymmInfo& info );
void DistSymmetricAnalysis
( const SymmElimTree& eTree, SymmInfo& info, 
  bool storeFactRecvIndices=true );

void ComputeFactRecvIndices
( const DistSymmNodeInfo& node,
  const DistSymmNodeInfo& childNode );

inline void SymmetricAnalysis
( const SymmElimTree& eTree, SymmInfo& info, bool storeFactRecvIndices )
{
#ifndef RELEASE
    PushCallStack("SymmetricAnalysis");
#endif
    LocalSymmetricAnalysis( eTree, info );
    DistSymmetricAnalysis( eTree, info, storeFactRecvIndices );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

#endif /* CLIQUE_SYMM_ANALYSIS_HPP */
