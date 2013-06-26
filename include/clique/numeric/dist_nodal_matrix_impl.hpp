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
DistNodalMatrix<F>::DistNodalMatrix()
: height_(0), width_(0)
{ }

template<typename F>
inline
DistNodalMatrix<F>::DistNodalMatrix
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::DistNodalMatrix");
#endif
    Pull( inverseMap, info, X );
}

template<typename F>
inline void
DistNodalMatrix<F>::Pull
( const DistMap& inverseMap, const DistSymmInfo& info,
  const DistMultiVec<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::Pull");
#endif
    throw std::logic_error("This must be rewritten");
}

template<typename F>
inline void
DistNodalMatrix<F>::Push
( const DistMap& inverseMap, const DistSymmInfo& info,
        DistMultiVec<F>& X ) const
{
#ifndef RELEASE
    CallStackEntry cse("DistNodalMatrix::Push");
#endif
    throw std::logic_error("This must be rewritten");
}

template<typename F>
inline int
DistNodalMatrix<F>::Height() const
{ return height_; }

template<typename F>
inline int
DistNodalMatrix<F>::Width() const
{ return width_; }

} // namespace cliq
