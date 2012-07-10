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
#ifndef CLIQUE_DIST_SEPARATOR_TREE_HPP
#define CLIQUE_DIST_SEPARATOR_TREE_HPP 1

namespace cliq {

struct LocalSepOrLeaf
{
    int localParent; // -1 if local root
    std::vector<int> indices;
};

struct DistSeparator
{
    mpi::Comm comm;
    std::vector<int> indices;
};

struct DistSeparatorTree
{
    // Full local binary tree
    std::vector<LocalSepOrLeaf> localSepsAndLeaves;
    // One path through top of binary tree (does not include local separator)
    std::vector<DistSeparator> distSeps;
};

} // namespace cliq

#endif /* CLIQUE_DIST_SEPARATOR_TREE_HPP */
