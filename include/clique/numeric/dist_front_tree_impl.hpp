/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename F>
inline
DistFrontTree<F>::DistFrontTree()
{ }

template<typename F>
inline
DistFrontTree<F>::DistFrontTree
( const DistSparseMatrix<F>& A, 
  const DistMap& map,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info )
: frontType(STRUCT_SYMM_2D)
{
#ifndef RELEASE
    CallStackEntry entry("DistFrontTree::DistFrontTree");
    if( A.LocalHeight() != map.NumLocalSources() )
        throw std::logic_error("Local mapping was not the right size");
#endif

    // TODO: Extend the DistSymmFrontTree version to handle the columns as well
    throw std::logic_error("This routine is not yet written");
}

} // namespace cliq
