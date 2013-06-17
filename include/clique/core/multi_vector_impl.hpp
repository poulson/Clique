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
MakeZeros( MultiVector<T>& X )
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
MakeUniform( MultiVector<T>& X )
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
Norms( const MultiVector<F>& X, std::vector<BASE(F)>& norms )
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

template<typename T>
inline void
Axpy( T alpha, const MultiVector<T>& X, MultiVector<T>& Y )
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
MultiVector<T>::MultiVector()
{ }

template<typename T>
inline
MultiVector<T>::MultiVector( int height, int width )
: multiVec_(height,width)
{ }

template<typename T>
inline 
MultiVector<T>::~MultiVector()
{ }

template<typename T>
inline int 
MultiVector<T>::Height() const
{ return multiVec_.Height(); }

template<typename T>
inline int
MultiVector<T>::Width() const
{ return multiVec_.Width(); }

template<typename T>
inline T
MultiVector<T>::Get( int row, int col ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("MultiVector::Get");
#endif
    return multiVec_.Get(row,col);
}

template<typename T>
inline void
MultiVector<T>::Set( int row, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("MultiVector::Set");
#endif
    multiVec_.Set(row,col,value);
}

template<typename T>
inline void
MultiVector<T>::Update( int row, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("MultiVector::Update");
#endif
    multiVec_.Update(row,col,value);
}

template<typename T>
inline void
MultiVector<T>::Empty()
{
    multiVec_.Empty();
}

template<typename T>
inline void
MultiVector<T>::ResizeTo( int height, int width )
{
    multiVec_.ResizeTo( height, width );
}

template<typename T>
const MultiVector<T>& 
MultiVector<T>::operator=( const Vector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("MultiVector::operator=");
#endif
    multiVec_ = x.vec_;
    return *this;
}

template<typename T>
const MultiVector<T>& 
MultiVector<T>::operator=( const MultiVector<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("MultiVector::operator=");
#endif
    multiVec_ = X.multiVec_;
    return *this;
}

} // namespace cliq
