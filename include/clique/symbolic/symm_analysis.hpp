/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

void SymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, 
  bool storeFactRecvIndices=true );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

void LocalSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info );
void DistSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, 
  bool storeFactRecvIndices=true );

void ComputeFactRecvIndices
( const DistSymmNodeInfo& node,
  const DistSymmNodeInfo& childNode );

inline void SymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool storeFactRecvIndices )
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
