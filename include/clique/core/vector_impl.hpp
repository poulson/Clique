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

namespace cliq {

template<typename T>
inline void 
MakeZeros( Vector<T>& x )
{
#ifndef RELEASE
    PushCallStack("MakeZeros");
#endif
    const int height = x.Height();
    for( int i=0; i<height; ++i )
        x.Set( i, T(0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
MakeUniform( Vector<T>& x )
{
#ifndef RELEASE
    PushCallStack("MakeUniform");
#endif
    const int height = x.Height();
    for( int i=0; i<height; ++i )
        x.Set( i, elem::SampleUnitBall<T>() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline typename Base<F>::type 
Norm( const Vector<F>& x )
{
#ifndef RELEASE
    PushCallStack("Norm");
#endif
    typedef typename Base<F>::type R;
    const int height = x.Height();

    R scale = 0;
    R scaledSquare = 1;
    for( int i=0; i<height; ++i )
    {
        const R alphaAbs = Abs(x.Get(i));
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
    const R norm = scale*Sqrt(scaledSquare);
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename T>
inline void
Axpy( T alpha, const Vector<T>& x, Vector<T>& y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
    if( x.Height() != y.Height() )
        throw std::logic_error("x and y must be the same height");
#endif
    const int height = x.Height(); 
    for( int i=0; i<height; ++i )
        y.Update( i, alpha*x.Get(i) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
Vector<T>::Vector()
{ }

template<typename T>
inline 
Vector<T>::Vector( int height )
: vec_(height,1)
{ }

template<typename T>
inline 
Vector<T>::~Vector()
{ }

template<typename T>
inline int 
Vector<T>::Height() const
{ return vec_.Height(); }

template<typename T>
inline T
Vector<T>::Get( int row ) const
{ 
#ifndef RELEASE 
    PushCallStack("Vector::Get");
#endif
    const T value = vec_.Get(row,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline void
Vector<T>::Set( int row, T value )
{
#ifndef RELEASE
    PushCallStack("Vector::Set");
#endif
    vec_.Set(row,0,value);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Vector<T>::Update( int row, T value )
{
#ifndef RELEASE
    PushCallStack("Vector::Update");
#endif
    vec_.Update(row,0,value);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Vector<T>::Empty()
{
    vec_.Empty();
}

template<typename T>
inline void
Vector<T>::ResizeTo( int height )
{
    vec_.ResizeTo( height, 1 );
}

template<typename T>
const Vector<T>& 
Vector<T>::operator=( const Vector<T>& x )
{
#ifndef RELEASE
    PushCallStack("Vector::operator=");
#endif
    vec_ = x.vec_;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace cliq
