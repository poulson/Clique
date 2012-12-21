/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

// For handling a vector distributed in a [VC,* ] manner over each node
// of the elimination tree
//
// TODO: Generalize to a set of vectors
//
template<typename F>
class DistNodalVector
{
public:
    Matrix<F> localVec;

    void Pull
    ( const DistMap& inverseMap, const DistSymmInfo& info,
      const DistVector<F>& x );
    void Push
    ( const DistMap& inverseMap, const DistSymmInfo& info,
            DistVector<F>& x ) const;

    DistNodalVector();
    DistNodalVector
    ( const DistMap& inverseMap, const DistSymmInfo& info, 
      const DistVector<F>& x );
};

} // namespace cliq
