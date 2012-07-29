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
MakeZeros( DistVector<T>& x )
{
#ifndef RELEASE
    PushCallStack("MakeZeros");
#endif
    const int localHeight = x.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        x.SetLocal( iLocal, (T)0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void 
MakeUniform( DistVector<T>& x )
{
#ifndef RELEASE
    PushCallStack("MakeUniform");
#endif
    const int localHeight = x.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        x.SetLocal( iLocal, elem::SampleUnitBall<T>() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline typename Base<F>::type 
Norm( const DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("Norm");
#endif
    typedef typename Base<F>::type R;
    const int localHeight = x.LocalHeight();
    mpi::Comm comm = x.Comm();

    R localScale = 0;
    R localScaledSquare = 1;
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        const R alphaAbs = Abs(x.GetLocal(iLocal));
        if( alphaAbs != 0 )
        {
            if( alphaAbs <= localScale )
            {
                const R relScale = alphaAbs/localScale;
                localScaledSquare += relScale*relScale;
            }
            else
            {
                const R relScale = localScale/alphaAbs;
                localScaledSquare = localScaledSquare*relScale*relScale + 1;
                localScale = alphaAbs;
            }
        }
    }

    // Find the maximum relative scale
    R scale;
    mpi::AllReduce( &localScale, &scale, 1, mpi::MAX, comm );

    R norm = 0;
    if( scale != 0 )
    {
        // Equilibrate our local scaled sum to the maximum scale
        R relScale = localScale/scale;
        localScaledSquare *= relScale*relScale;

        // The scaled square is now simply the sum of the local contributions
        R scaledSquare;
        mpi::AllReduce( &localScaledSquare, &scaledSquare, 1, mpi::SUM, comm );

        norm = scale*Sqrt(scaledSquare);
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename T>
inline void
Axpy( T alpha, const DistVector<T>& x, DistVector<T>& y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
    if( !mpi::CongruentComms( x.Comm(), y.Comm() ) )
        throw std::logic_error("x and y must have congruent communicators");
    if( x.Height() != y.Height() )
        throw std::logic_error("x and y must be the same height");
#endif
    const int localHeight = x.LocalHeight(); 
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        y.UpdateLocal( iLocal, alpha*x.GetLocal(iLocal) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
DistVector<T>::DistVector()
: height_(0), comm_(mpi::COMM_WORLD), blocksize_(0), firstLocalRow_(0)
{ }

template<typename T>
inline 
DistVector<T>::DistVector( mpi::Comm comm )
: height_(0), blocksize_(0), firstLocalRow_(0)
{ 
    SetComm( comm );
}

template<typename T>
inline 
DistVector<T>::DistVector( int height, mpi::Comm comm )
: height_(height)
{ 
    SetComm( comm );
}

template<typename T>
inline 
DistVector<T>::~DistVector()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommFree( comm_ );
}

template<typename T>
inline int 
DistVector<T>::Height() const
{ return height_; }

template<typename T>
inline void
DistVector<T>::SetComm( mpi::Comm comm )
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommDup( comm, comm_ );
    else
        comm_ = comm;

    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    blocksize_ = height_/commSize;
    firstLocalRow_ = blocksize_*commRank;
    const int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    localVec_.ResizeTo( localHeight, 1 );
}

template<typename T>
inline mpi::Comm 
DistVector<T>::Comm() const
{ return comm_; }

template<typename T>
inline int
DistVector<T>::Blocksize() const
{ return blocksize_; }

template<typename T>
inline int
DistVector<T>::FirstLocalRow() const
{ return firstLocalRow_; }

template<typename T>
inline int
DistVector<T>::LocalHeight() const
{ return localVec_.Height(); }

template<typename T>
inline T
DistVector<T>::GetLocal( int localRow ) const
{ 
#ifndef RELEASE 
    PushCallStack("DistVector::GetLocal");
#endif
    const T value = localVec_.Get(localRow,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline void
DistVector<T>::SetLocal( int localRow, T value )
{
#ifndef RELEASE
    PushCallStack("DistVector::SetLocal");
#endif
    localVec_.Set(localRow,0,value);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistVector<T>::UpdateLocal( int localRow, T value )
{
#ifndef RELEASE
    PushCallStack("DistVector::UpdateLocal");
#endif
    localVec_.Update(localRow,0,value);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DistVector<T>::Empty()
{
    height_ = 0;
    blocksize_ = 0;
    firstLocalRow_ = 0;
    localVec_.Empty();
}

template<typename T>
inline void
DistVector<T>::ResizeTo( int height )
{
    const int commRank = mpi::CommRank( comm_ );
    const int commSize = mpi::CommSize( comm_ );
    height_ = height;
    blocksize_ = height/commSize;
    firstLocalRow_ = commRank*blocksize_;
    const int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    localVec_.ResizeTo( localHeight, 1 );
}

template<typename T>
const DistVector<T>& 
DistVector<T>::operator=( const DistVector<T>& x )
{
#ifndef RELEASE
    PushCallStack("DistVector::operator=");
#endif
    height_ = x.height_;
    SetComm( x.comm_ );
    localVec_ = x.localVec_;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace cliq
