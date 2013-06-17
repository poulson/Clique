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
MakeZeros( DistMultiVector<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("MakeZeros");
#endif
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            X.SetLocal( iLocal, j, T(0) );
}

template<typename T>
inline void 
MakeUniform( DistMultiVector<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("MakeUniform");
#endif
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            X.SetLocal( iLocal, j, elem::SampleUnitBall<T>() );
}

template<typename F>
void Norms( const DistMultiVector<F>& X, std::vector<BASE(F)>& norms )
{
#ifndef RELEASE
    CallStackEntry entry("Norms");
#endif
    typedef BASE(F) R;
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    mpi::Comm comm = X.Comm();

    norms.resize( width );
    std::vector<R> localScales( width ),
                   localScaledSquares( width );
    for( int j=0; j<width; ++j )
    {
        R localScale = 0;
        R localScaledSquare = 1;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const R alphaAbs = Abs(X.GetLocal(iLocal,j));
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

        localScales[j] = localScale;
        localScaledSquares[j] = localScaledSquare;
    }

    // Find the maximum relative scales
    std::vector<R> scales( width );
    mpi::AllReduce( &localScales[0], &scales[0], width, mpi::MAX, comm );

    // Equilibrate the local scaled sums
    for( int j=0; j<width; ++j )
    {
        const R scale = scales[j];
        if( scale != 0 )
        {
            // Equilibrate our local scaled sum to the maximum scale
            R relScale = localScales[j]/scale;
            localScaledSquares[j] *= relScale*relScale;
        }
        else
            localScaledSquares[j] = 0;
    }

    // Combine the local contributions
    std::vector<R> scaledSquares( width );
    mpi::AllReduce
    ( &localScaledSquares[0], &scaledSquares[0], width, mpi::SUM, comm );
    for( int j=0; j<width; ++j )
        norms[j] = scales[j]*Sqrt(scaledSquares[j]);
}

template<typename T>
inline void
Axpy( T alpha, const DistMultiVector<T>& X, DistMultiVector<T>& Y )
{
#ifndef RELEASE
    CallStackEntry entry("Axpy");
    if( !mpi::CongruentComms( X.Comm(), Y.Comm() ) )
        throw std::logic_error("X and Y must have congruent communicators");
    if( X.Height() != Y.Height() )
        throw std::logic_error("X and Y must be the same height");
    if( X.Width() != Y.Width() )
        throw std::logic_error("X and Y must be the same width");
#endif
    const int localHeight = X.LocalHeight(); 
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            Y.UpdateLocal( iLocal, j, alpha*X.GetLocal(iLocal,j) );
}

template<typename T>
inline 
DistMultiVector<T>::DistMultiVector()
: height_(0), comm_(mpi::COMM_WORLD), blocksize_(0), firstLocalRow_(0)
{ }

template<typename T>
inline 
DistMultiVector<T>::DistMultiVector( mpi::Comm comm )
: height_(0), width_(0), comm_(mpi::COMM_WORLD), 
  blocksize_(0), firstLocalRow_(0)
{ 
    SetComm( comm );
}

template<typename T>
inline 
DistMultiVector<T>::DistMultiVector( int height, int width, mpi::Comm comm )
: height_(height), width_(width), comm_(mpi::COMM_WORLD)
{ 
    SetComm( comm );
}

template<typename T>
inline 
DistMultiVector<T>::~DistMultiVector()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommFree( comm_ );
}

template<typename T>
inline int 
DistMultiVector<T>::Height() const
{ return height_; }

template<typename T>
inline int
DistMultiVector<T>::Width() const
{ return multiVec_.Width(); }

template<typename T>
inline void
DistMultiVector<T>::SetComm( mpi::Comm comm )
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
    multiVec_.ResizeTo( localHeight, width_ );
}

template<typename T>
inline mpi::Comm 
DistMultiVector<T>::Comm() const
{ return comm_; }

template<typename T>
inline int
DistMultiVector<T>::Blocksize() const
{ return blocksize_; }

template<typename T>
inline int
DistMultiVector<T>::FirstLocalRow() const
{ return firstLocalRow_; }

template<typename T>
inline int
DistMultiVector<T>::LocalHeight() const
{ return multiVec_.Height(); }

template<typename T>
inline T
DistMultiVector<T>::GetLocal( int localRow, int col ) const
{ 
#ifndef RELEASE 
    CallStackEntry entry("DistMultiVector::GetLocal");
#endif
    return multiVec_.Get(localRow,col);
}

template<typename T>
inline void
DistMultiVector<T>::SetLocal( int localRow, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("DistMultiVector::SetLocal");
#endif
    multiVec_.Set(localRow,col,value);
}

template<typename T>
inline void
DistMultiVector<T>::UpdateLocal( int localRow, int col, T value )
{
#ifndef RELEASE
    CallStackEntry entry("DistMultiVector::UpdateLocal");
#endif
    multiVec_.Update(localRow,col,value);
}

template<typename T>
inline void
DistMultiVector<T>::Empty()
{
    height_ = 0;
    width_ = 0;
    blocksize_ = 0;
    firstLocalRow_ = 0;
    multiVec_.Empty();
}

template<typename T>
inline void
DistMultiVector<T>::ResizeTo( int height, int width )
{
    const int commRank = mpi::CommRank( comm_ );
    const int commSize = mpi::CommSize( comm_ );
    height_ = height;
    width_ = width;
    blocksize_ = height/commSize;
    firstLocalRow_ = commRank*blocksize_;
    const int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    multiVec_.ResizeTo( localHeight, width );
}

template<typename T>
const DistMultiVector<T>& 
DistMultiVector<T>::operator=( const DistVector<T>& x )
{
#ifndef RELEASE
    CallStackEntry entry("DistMultiVector::operator=");
#endif
    height_ = x.height_;
    width_ = 1;
    SetComm( x.comm_ );
    multiVec_ = x.localVec_;
    return *this;
}

template<typename T>
const DistMultiVector<T>& 
DistMultiVector<T>::operator=( const DistMultiVector<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DistMultiVector::operator=");
#endif
    height_ = X.height_;
    width_ = X.width_;
    SetComm( X.comm_ );
    multiVec_ = X.multiVec_;
    return *this;
}

} // namespace cliq
