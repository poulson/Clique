/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_SYMBOLIC_SYMMANALYSIS_HPP
#define CLIQ_SYMBOLIC_SYMMANALYSIS_HPP

namespace cliq {

void SymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, 
  bool storeFactRecvInd=true );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

void LocalSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info );
void DistSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, 
  bool storeFactRecvInd=true );

inline void SymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool storeFactRecvInd )
{
#ifndef RELEASE
    CallStackEntry entry("SymmetricAnalysis");
#endif
    LocalSymmetricAnalysis( eTree, info );
    DistSymmetricAnalysis( eTree, info, storeFactRecvInd );
}

} // namespace cliq

#endif // ifndef CLIQ_SYMBOLIC_SYMMANALYSIS_HPP
