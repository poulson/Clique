/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_SYMBOLIC_DISTSYMMINFO_HPP
#define CLIQ_SYMBOLIC_DISTSYMMINFO_HPP

namespace cliq {

struct SymmNodeInfo
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
    bool onLeft;
    int myOffset;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelInd;
    // (maps from the child update indices to our frontal indices).
    std::vector<int> leftRelInd, rightRelInd;
};

struct FactorCommMeta
{
    std::vector<int> numChildSendInd;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::vector<int> > childRecvInd;

    void EmptyChildRecvIndices() const
    { std::vector<std::vector<int> >().swap(childRecvInd); }

    void Empty()
    {
        std::vector<int>().swap(numChildSendInd);
        EmptyChildRecvIndices();
    }
};

struct MultiVecCommMeta
{
    int localSize;
    std::vector<int> numChildSendInd;
    std::vector<std::vector<int> > childRecvInd;

    void Empty()
    {
        std::vector<int>().swap( numChildSendInd );
        std::vector<std::vector<int> >().swap( childRecvInd );
    }
};

struct MatrixCommMeta
{
    std::vector<int> numChildSendInd;
    std::vector<std::vector<int> > childRecvInd;

    void Empty()
    {
        std::vector<int>().swap( numChildSendInd );
        std::vector<std::vector<int> >().swap( childRecvInd );
    }
};

struct DistSymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, offset;
    std::vector<int> origLowerStruct;
    bool onLeft;
    mpi::Comm comm;

    //
    // The following is computed during analysis
    //
    Grid* grid;
    int myOffset, leftSize, rightSize;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelInd;

    // The relative indices of our child
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftRelInd, rightRelInd;

    FactorCommMeta factorMeta;
    MultiVecCommMeta multiVecMeta;
};

struct DistSymmInfo
{
    std::vector<SymmNodeInfo> localNodes;
    std::vector<DistSymmNodeInfo> distNodes;
    ~DistSymmInfo();
};

// Utilities
void ComputeFactRecvInd
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode );
void GetChildGridDims
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode,
  int* childGridDims );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline
DistSymmInfo::~DistSymmInfo()
{
    const int numDist = distNodes.size();
    for( int s=0; s<numDist; ++s )
    {
        delete distNodes[s].grid;
        mpi::CommFree( distNodes[s].comm );
    }
}

} // namespace cliq

#endif // ifndef CLIQ_SYMBOLIC_DISTSYMMINFO_HPP
