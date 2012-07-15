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
#ifndef CLIQUE_HPP
#define CLIQUE_HPP 1

#include "clique/core/environment.hpp"

#include "clique/symbolic/dist_graph.hpp"
#include "clique/symbolic/graph.hpp"
#include "clique/symbolic/dist_separator_tree.hpp"
#include "clique/symbolic/dist_symm_elim_tree.hpp"
#include "clique/symbolic/dist_symm_info.hpp"
#include "clique/symbolic/symm_analysis.hpp"
#include "clique/symbolic/nested_dissection.hpp"

// These two classes are interdependent, so declare both, then define both
#include "clique/numeric/dist_vector.hpp"
#include "clique/numeric/dist_nodal_vector.hpp"
#include "clique/numeric/dist_vector_main.hpp"
#include "clique/numeric/dist_nodal_vector_main.hpp"

#include "clique/numeric/dist_sparse_matrix.hpp"
#include "clique/numeric/sparse_matrix.hpp"
#include "clique/numeric/multiply.hpp"

#include "clique/numeric/dist_symm_front_tree.hpp"
#include "clique/numeric/set_solve_mode.hpp"

#include "clique/numeric/ldl.hpp"
#include "clique/numeric/ldl_solve.hpp"
#include "clique/numeric/lower_solve.hpp"
#include "clique/numeric/diagonal_solve.hpp"
#include "clique/numeric/lower_multiply.hpp"

#include "clique/numeric/block_lower_solve.hpp"
#include "clique/numeric/block_ldl.hpp"
#include "clique/numeric/block_ldl_solve.hpp"

#endif /* CLIQUE_HPP */
