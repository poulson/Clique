/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/
#ifndef CLIQUE_HPP
#define CLIQUE_HPP 1

#include "elemental.hpp"
#include <algorithm>
#include <map>
#include <set>

#include "clique/config.h"

//
// The core of the library
//

#include "clique/core/environment.hpp"

// Data-structure declarations
namespace cliq { 
    class Graph; 
    class DistGraph; 
    template<typename T> class SparseMatrix;
    template<typename T> class DistSparseMatrix;
    template<typename T> class Vector;
    template<typename T> class DistVector;
} 
#include "clique/core/graph_decl.hpp"
#include "clique/core/dist_graph_decl.hpp"
//TODO: include "clique/core/map_decl.hpp"
#include "clique/core/dist_map_decl.hpp"
#include "clique/core/sparse_matrix_decl.hpp"
#include "clique/core/dist_sparse_matrix_decl.hpp"
#include "clique/core/vector_decl.hpp"
#include "clique/core/dist_vector_decl.hpp"
#include "clique/core/multi_vector_decl.hpp"
#include "clique/core/dist_multi_vector_decl.hpp"

// Data-structure implementations
#include "clique/core/graph_impl.hpp"
#include "clique/core/dist_graph_impl.hpp"
//TODO: include "clique/core/map_impl.hpp"
#include "clique/core/dist_map_impl.hpp"
#include "clique/core/sparse_matrix_impl.hpp"
#include "clique/core/dist_sparse_matrix_impl.hpp"
#include "clique/core/vector_impl.hpp"
#include "clique/core/dist_vector_impl.hpp"
#include "clique/core/multi_vector_impl.hpp"
#include "clique/core/dist_multi_vector_impl.hpp"

//
// Symbolic computation
//

#include "clique/symbolic/dist_separator_tree.hpp"
#include "clique/symbolic/dist_symm_elim_tree.hpp"
#include "clique/symbolic/dist_symm_info.hpp"
#include "clique/symbolic/symm_analysis.hpp"
#include "clique/symbolic/nested_dissection.hpp"
#include "clique/symbolic/natural_nested_dissection.hpp"

//
// Numerical computation
//

#include "clique/numeric/multiply.hpp"

#include "clique/numeric/dist_nodal_vector_decl.hpp"
#include "clique/numeric/dist_nodal_vector_impl.hpp"

#include "clique/numeric/dist_nodal_multi_vector_decl.hpp"
#include "clique/numeric/dist_nodal_multi_vector_impl.hpp"

#include "clique/numeric/dist_front_tree_decl.hpp"
#include "clique/numeric/dist_front_tree_impl.hpp"
#include "clique/numeric/dist_symm_front_tree_decl.hpp"
#include "clique/numeric/dist_symm_front_tree_impl.hpp"
#include "clique/numeric/change_front_type.hpp"

#include "clique/numeric/ldl.hpp"
#include "clique/numeric/lower_solve.hpp"
#include "clique/numeric/diagonal_solve.hpp"
#include "clique/numeric/solve.hpp"
#include "clique/numeric/lower_multiply.hpp"

#endif // CLIQUE_HPP
