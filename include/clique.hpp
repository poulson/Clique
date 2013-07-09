/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef CLIQUE_HPP
#define CLIQUE_HPP

#include "elemental.hpp"
#include <algorithm>
#include <map>
#include <set>

#include "clique/config.h"

//
// The core of the library
//

#include "clique/core/environment.hpp"

// Data-structure declarations (are these needed?!?)
namespace cliq { 
    class Graph; 
    class DistGraph; 
    template<typename T> class SparseMatrix;
    template<typename T> class DistSparseMatrix;
    template<typename T> class DistNodalMultiVec;
    template<typename T> class DistNodalMatrix;
} 
#include "clique/core/graph_decl.hpp"
#include "clique/core/dist_graph_decl.hpp"
//TODO: include "clique/core/map_decl.hpp"
#include "clique/core/dist_map_decl.hpp"
#include "clique/core/sparse_matrix_decl.hpp"
#include "clique/core/dist_sparse_matrix_decl.hpp"
#include "clique/core/multi_vec_decl.hpp"
#include "clique/core/dist_multi_vec_decl.hpp"

// Data-structure implementations
#include "clique/core/graph_impl.hpp"
#include "clique/core/dist_graph_impl.hpp"
//TODO: include "clique/core/map_impl.hpp"
#include "clique/core/dist_map_impl.hpp"
#include "clique/core/sparse_matrix_impl.hpp"
#include "clique/core/dist_sparse_matrix_impl.hpp"
#include "clique/core/multi_vec_impl.hpp"
#include "clique/core/dist_multi_vec_impl.hpp"

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

#include "clique/numeric/dist_nodal_multi_vec_decl.hpp"
#include "clique/numeric/dist_nodal_matrix_decl.hpp"
#include "clique/numeric/dist_nodal_multi_vec_impl.hpp"
#include "clique/numeric/dist_nodal_matrix_impl.hpp"

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

//
// Input/Output
//

#include "clique/io/Print.hpp"
#include "clique/io/Display.hpp"

#endif // ifndef CLIQUE_HPP
