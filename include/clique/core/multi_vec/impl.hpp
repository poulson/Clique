/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_CORE_MULTIVEC_IMPL_HPP
#define CLIQ_CORE_MULTIVEC_IMPL_HPP

namespace cliq {

template<typename T>
inline void 
MakeZeros( MultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MakeZeros"))
    const int height = X.Height();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            X.Set( i, j, T(0) );
}

template<typename T>
inline void 
MakeUniform( MultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    const int height = X.Height();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            X.Set( i, j, elem::SampleBall<T>() );
}

template<typename F>
inline void
Norms( const MultiVec<F>& X, std::vector<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("Norms"))
    typedef BASE(F) R;
    const int height = X.Height();
    const int width = X.Width();

    norms.resize( width );
    for( int j=0; j<width; ++j )
    {
        R scale = 0;
        R scaledSquare = 1;
        for( int i=0; i<height; ++i )
        {
            const R alphaAbs = Abs(X.Get(i,j));
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= scale )
                {
                    const R relScale = alphaAbs/scale;
                    scaledSquare += relScale*relScale;
                }
                else
                {
                    const R relScale = scale/alphaAbs;
                    scaledSquare = scaledSquare*relScale*relScale + 1;
                    scale = alphaAbs;
                }
            }
        }
        norms[j] = scale*Sqrt(scaledSquare);
    }
}

template<typename F>
inline BASE(F)
Norm( const MultiVec<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Norm"))
    if( x.Width() != 1 )
        LogicError("Norms only applies with one column");
    typedef BASE(F) R;
    std::vector<R> norms;
    Norms( x, norms );
    return norms[0];
}

template<typename T>
inline void
Axpy( T alpha, const MultiVec<T>& X, MultiVec<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Axpy");
        if( X.Height() != Y.Height() )
            LogicError("X and Y must be the same height");
        if( X.Width() != Y.Width() )
            LogicError("X and Y must be the same width");
    )
    const int height = X.Height(); 
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            Y.Update( i, j, alpha*X.Get(i,j) );
}

template<typename T>
inline 
MultiVec<T>::MultiVec()
{ }

template<typename T>
inline
MultiVec<T>::MultiVec( int height, int width )
: multiVec_(height,width)
{ }

template<typename T>
inline 
MultiVec<T>::~MultiVec()
{ }

template<typename T>
inline int 
MultiVec<T>::Height() const
{ return multiVec_.Height(); }

template<typename T>
inline int
MultiVec<T>::Width() const
{ return multiVec_.Width(); }

template<typename T>
inline T
MultiVec<T>::Get( int row, int col ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("MultiVec::Get"))
    return multiVec_.Get(row,col);
}

template<typename T>
inline void
MultiVec<T>::Set( int row, int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("MultiVec::Set"))
    multiVec_.Set(row,col,value);
}

template<typename T>
inline void
MultiVec<T>::Update( int row, int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("MultiVec::Update"))
    multiVec_.Update(row,col,value);
}

template<typename T>
inline void
MultiVec<T>::Empty()
{ multiVec_.Empty(); }

template<typename T>
inline void
MultiVec<T>::ResizeTo( int height, int width )
{ multiVec_.ResizeTo( height, width ); }

template<typename T>
const MultiVec<T>& 
MultiVec<T>::operator=( const MultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiVec::operator="))
    multiVec_ = X.multiVec_;
    return *this;
}

} // namespace cliq

#endif // ifndef CLIQ_CORE_MULTIVEC_IMPL_HPP
