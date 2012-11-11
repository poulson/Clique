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

namespace cliq {

enum SymmFrontType
{
  SYMM_1D,
  SYMM_2D,
  LDL_1D,
  LDL_2D,
  LDL_SELINV_1D,
  LDL_SELINV_2D,
  BLOCK_LDL_2D
};

inline bool
FrontsAre1d( SymmFrontType frontType )
{
    if( frontType == SYMM_1D ||
        frontType == LDL_1D || 
        frontType == LDL_SELINV_1D )
        return true;
    else
        return false;
}

// Only keep track of the left and bottom-right piece of the fronts
// (with the bottom-right piece stored in workspace) since only the left side
// needs to be kept after the factorization is complete.

template<typename F>
struct LocalSymmFront
{
    Matrix<F> frontL;
    mutable Matrix<F> work;
};

template<typename F>
struct DistSymmFront
{
    // The 'frontType' member variable of the parent 'DistSymmFrontTree' 
    // determines which of the following fronts is active.
    //
    // Split each front into a left and right piece such that the right piece
    // is not needed after the factorization (and can be freed).

    // TODO: Think about the fact that almost everything is now mutable...

    mutable DistMatrix<F,VC,STAR> front1dL;
    mutable DistMatrix<F,VC,STAR> work1d;

    mutable DistMatrix<F> front2dL;
    mutable DistMatrix<F> work2d;

    DistMatrix<F,VC,STAR> diag;
};

template<typename F>
struct DistSymmFrontTree
{
    bool isHermitian;
    SymmFrontType frontType;
    std::vector<LocalSymmFront<F> > localFronts;
    std::vector<DistSymmFront<F> > distFronts;

    DistSymmFrontTree();

    DistSymmFrontTree
    ( Orientation orientation,
      const DistSparseMatrix<F>& A, 
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info );

    void Initialize
    ( Orientation orientation,
      const DistSparseMatrix<F>& A,
      const DistMap& reordering,
      const DistSeparatorTree& sepTree,
      const DistSymmInfo& info );

    void TopLeftMemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void BottomLeftMemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;

    void MemoryInfo
    ( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries,
      double& numGlobalEntries ) const;
};

} // namespace cliq
