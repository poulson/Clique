/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
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
    PushCallStack("DistFrontTree::DistFrontTree");
    if( A.LocalHeight() != map.NumLocalSources() )
        throw std::logic_error("Local mapping was not the right size");
#endif
    mpi::Comm comm = A.Comm();
    const DistGraph& graph = A.Graph();
    const int blocksize = A.Blocksize();
    const int commSize = mpi::CommSize( comm );
    const int numSources = graph.NumSources();
    const int numLocal = sepTree.localSepsAndLeaves.size();
    const int numDist = sepTree.distSeps.size();

    // TODO: Extend the DistSymmFrontTree version to handle the columns as well
    throw std::logic_error("This routine is not yet written");

#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
