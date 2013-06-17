/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

// For handling a vector distributed in a [VC,* ] manner over each node
// of the elimination tree
//
// TODO: Generalize to a set of vectors
//
template<typename F>
class DistNodalMultiVector
{
public:
    Matrix<F> multiVec;

    void Pull
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVector<F>& X );
    void Push
    ( const DistMap& inverseMap, const DistSymmInfo& info,
            DistMultiVector<F>& X ) const;

    DistNodalMultiVector();
    DistNodalMultiVector
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistMultiVector<F>& X );
};

} // namespace cliq
