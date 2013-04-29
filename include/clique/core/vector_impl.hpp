/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T>
inline void 
MakeZeros( Vector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("MakeZeros");
#endif
    const int height = x.Height();
    for( int i=0; i<height; ++i )
        x.Set( i, T(0) );
}

template<typename T>
inline void 
MakeUniform( Vector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("MakeUniform");
#endif
    const int height = x.Height();
    for( int i=0; i<height; ++i )
        x.Set( i, elem::SampleUnitBall<T>() );
}

template<typename F>
inline BASE(F)
Norm( const Vector<F>& x )
{
#ifndef RELEASE
    CallStackEntry entry("Norm");
#endif
    typedef BASE(F) R;
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
    return scale*Sqrt(scaledSquare);
}

template<typename T>
inline void
Axpy( T alpha, const Vector<T>& x, Vector<T>& y )
{
#ifndef RELEASE
    CallStackEntry entry("Axpy");
    if( x.Height() != y.Height() )
        throw std::logic_error("x and y must be the same height");
#endif
    const int height = x.Height(); 
    for( int i=0; i<height; ++i )
        y.Update( i, alpha*x.Get(i) );
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
    CallStackEntry entry("Vector::Get");
#endif
    return vec_.Get(row,0);
}

template<typename T>
inline void
Vector<T>::Set( int row, T value )
{
#ifndef RELEASE
    CallStackEntry entry("Vector::Set");
#endif
    vec_.Set(row,0,value);
}

template<typename T>
inline void
Vector<T>::Update( int row, T value )
{
#ifndef RELEASE
    CallStackEntry entry("Vector::Update");
#endif
    vec_.Update(row,0,value);
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
    CallStackEntry entry("Vector::operator=");
#endif
    vec_ = x.vec_;
    return *this;
}

} // namespace cliq
