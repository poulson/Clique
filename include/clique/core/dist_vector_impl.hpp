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
MakeZeros( DistVector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("MakeZeros");
#endif
    const int localHeight = x.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        x.SetLocal( iLocal, T(0) );
}

template<typename T>
inline void 
MakeUniform( DistVector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("MakeUniform");
#endif
    const int localHeight = x.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        x.SetLocal( iLocal, elem::SampleUnitBall<T>() );
}

template<typename F>
inline BASE(F)
Norm( const DistVector<F>& x )
{
#ifndef RELEASE
    CallStackEntry entry("Norm");
#endif
    typedef BASE(F) R;
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
    return norm;
}

template<typename T>
inline void
Axpy( T alpha, const DistVector<T>& x, DistVector<T>& y )
{
#ifndef RELEASE
    CallStackEntry entry("Axpy");
    if( !mpi::CongruentComms( x.Comm(), y.Comm() ) )
        throw std::logic_error("x and y must have congruent communicators");
    if( x.Height() != y.Height() )
        throw std::logic_error("x and y must be the same height");
#endif
    const int localHeight = x.LocalHeight(); 
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        y.UpdateLocal( iLocal, alpha*x.GetLocal(iLocal) );
}

template<typename T>
inline 
DistVector<T>::DistVector()
: height_(0), comm_(mpi::COMM_WORLD), blocksize_(0), firstLocalRow_(0)
{ }

template<typename T>
inline 
DistVector<T>::DistVector( mpi::Comm comm )
: height_(0), comm_(mpi::COMM_WORLD), blocksize_(0), firstLocalRow_(0)
{ 
    SetComm( comm );
}

template<typename T>
inline 
DistVector<T>::DistVector( int height, mpi::Comm comm )
: height_(height), comm_(mpi::COMM_WORLD)
{ 
    SetComm( comm );
}

template<typename T>
inline 
DistVector<T>::DistVector( int height, T* buffer, mpi::Comm comm )
: height_(height)
{
    if( comm != mpi::COMM_WORLD )
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
    vec_.Attach( localHeight, 1, buffer, localHeight );
}

template<typename T>
inline 
DistVector<T>::DistVector( int height, const T* buffer, mpi::Comm comm )
: height_(height)
{
    if( comm != mpi::COMM_WORLD )
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
    vec_.LockedAttach( localHeight, 1, buffer, localHeight );
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
        mpi::CommFree( comm_ );

    if( comm != mpi::COMM_WORLD )
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
    vec_.ResizeTo( localHeight, 1 );
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
{ return vec_.Height(); }

template<typename T>
inline T
DistVector<T>::GetLocal( int localRow ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistVector::GetLocal");
#endif
    return vec_.Get(localRow,0);
}

template<typename T>
inline void
DistVector<T>::SetLocal( int localRow, T value )
{
#ifndef RELEASE
    CallStackEntry entry("DistVector::SetLocal");
#endif
    vec_.Set(localRow,0,value);
}

template<typename T>
inline void
DistVector<T>::UpdateLocal( int localRow, T value )
{
#ifndef RELEASE
    CallStackEntry entry("DistVector::UpdateLocal");
#endif
    vec_.Update(localRow,0,value);
}

template<typename T>
inline Matrix<T>&
DistVector<T>::Vector()
{ return vec_; }

template<typename T>
inline const Matrix<T>&
DistVector<T>::Vector() const
{ return vec_; }

template<typename T>
inline void
DistVector<T>::Empty()
{
    height_ = 0;
    blocksize_ = 0;
    firstLocalRow_ = 0;
    vec_.Empty();
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
    vec_.ResizeTo( localHeight, 1 );
}

template<typename T>
const DistVector<T>& 
DistVector<T>::operator=( const DistVector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("DistVector::operator=");
#endif
    height_ = x.height_;
    SetComm( x.comm_ );
    vec_ = x.vec_;
    return *this;
}

} // namespace cliq
