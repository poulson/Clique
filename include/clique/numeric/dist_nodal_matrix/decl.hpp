/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_DISTNODALMATRIX_DECL_HPP
#define CLIQ_NUMERIC_DISTNODALMATRIX_DECL_HPP

namespace cliq {

// For handling a matrix distributed in a [MC,MR] manner over each node
// of the elimination tree
template<typename F>
class DistNodalMatrix
{
public:
    std::vector<Matrix<F>> localNodes;
    std::vector<DistMatrix<F>> distNodes;

    DistNodalMatrix();
    DistNodalMatrix
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVec<F>& X );
    DistNodalMatrix( const DistNodalMultiVec<F>& X );

    const DistNodalMatrix<F>& operator=( const DistNodalMultiVec<F>& X );

    void Pull
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVec<F>& X );
    void Push
    ( const DistMap& inverseMap, const DistSymmInfo& info,
            DistMultiVec<F>& X ) const;

    int Height() const;
    int Width() const;

    mutable std::vector<MatrixCommMeta> commMetas;
    void ComputeCommMetas( const DistSymmInfo& info ) const;
private:
    int height_, width_;
};

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_DISTNODALMATRIX_DECL_HPP
