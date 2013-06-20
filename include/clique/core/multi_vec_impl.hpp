/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename T>
inline void 
MakeZeros( MultiVec<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("MakeZeros");
#endif
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
#ifndef RELEASE
    CallStackEntry entry("MakeUniform");
#endif
    const int height = X.Height();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            X.Set( i, j, elem::SampleUnitBall<T>() );
}

template<typename F>
inline void
Norms( const MultiVec<F>& X, std::vector<BASE(F)>& norms )
{
#ifndef RELEASE
    CallStackEntry entry("Norms");
#endif
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
#ifndef RELEASE
    CallStackEntry entry("Norm");
#endif
    if( x.Width() != 1 )
        throw std::logic_error("Norms only applies with one column");
    typedef BASE(F) R;
    std::vector<R> norms;
    Norms( x, norms );
    return norms[0];
}

template<typename T>
inline void
Axpy( T alpha, const MultiVec<T>& X, MultiVec<T>& Y )
{
#ifndef RELEASE
    CallStackEntry entry("Axpy");
    if( X.Height() != Y.Height() )
        throw std::logic_error("X and Y must be the same height");
    if( X.Width() != Y.Width() )
        throw std::logic_error("X and Y must be the same width");
#endif
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
#ifndef RELEASE 
    CallStackEntry entry("MultiVec::Get");
#endif
    return multiVec_.Get(row,col);
}

template<typename T>
inline void
MultiVec<T>::Set( int row, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("MultiVec::Set");
#endif
    multiVec_.Set(row,col,value);
}

template<typename T>
inline void
MultiVec<T>::Update( int row, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("MultiVec::Update");
#endif
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
#ifndef RELEASE
    CallStackEntry entry("MultiVec::operator=");
#endif
    multiVec_ = X.multiVec_;
    return *this;
}

} // namespace cliq
