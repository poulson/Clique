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
#ifndef CLIQUE_DIST_VECTOR_MAIN_HPP
#define CLIQUE_DIST_VECTOR_MAIN_HPP 1

namespace cliq {

template<typename F>
void MakeZeros( DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("MakeZeros");
#endif
    const int localHeight = x.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        x.SetLocal( iLocal, (F)0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void MakeUniform( DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("MakeUniform");
#endif
    const int localHeight = x.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        x.SetLocal( iLocal, elem::SampleUnitBall<F>() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline 
DistVector<F>::DistVector()
: height_(0), comm_(mpi::COMM_WORLD), blocksize_(0), firstLocalRow_(0)
{ }

template<typename F>
inline 
DistVector<F>::DistVector( mpi::Comm comm )
: height_(0), blocksize_(0), firstLocalRow_(0)
{ 
    SetComm( comm );
}

template<typename F>
inline 
DistVector<F>::DistVector( int height, mpi::Comm comm )
: height_(height)
{ 
    SetComm( comm );
}

template<typename F>
inline 
DistVector<F>::~DistVector()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::CommFree( comm_ );
}

template<typename F>
inline int 
DistVector<F>::Height() const
{ return height_; }

template<typename F>
inline void
DistVector<F>::SetComm( mpi::Comm comm )
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
    values_.resize( localHeight );
}

template<typename F>
inline mpi::Comm 
DistVector<F>::Comm() const
{ return comm_; }

template<typename F>
inline int
DistVector<F>::Blocksize() const
{ return blocksize_; }

template<typename F>
inline int
DistVector<F>::FirstLocalRow() const
{ return firstLocalRow_; }

template<typename F>
inline int
DistVector<F>::LocalHeight() const
{ return values_.size(); }

template<typename F>
inline F
DistVector<F>::GetLocal( int localRow ) const
{ 
#ifndef RELEASE 
    PushCallStack("DistVector::GetLocal");
    if( localRow < 0 || localRow >= (int)values_.size() )
        throw std::logic_error("Local row out of bounds");
    PopCallStack();
#endif
    return values_[localRow];
}

template<typename F>
inline void
DistVector<F>::SetLocal( int localRow, F value )
{
#ifndef RELEASE
    PushCallStack("DistVector::SetLocal");
    if( localRow < 0 || localRow >= (int)values_.size() )
        throw std::logic_error("Local row out of bounds");
#endif
    values_[localRow] = value;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistVector<F>::UpdateLocal( int localRow, F value )
{
#ifndef RELEASE
    PushCallStack("DistVector::UpdateLocal");
    if( localRow < 0 || localRow >= (int)values_.size() )
        throw std::logic_error("Local row out of bounds");
#endif
    values_[localRow] += value;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
DistVector<F>::Empty()
{
    height_ = 0;
    blocksize_ = 0;
    firstLocalRow_ = 0;
    values_.clear();
}

template<typename F>
inline void
DistVector<F>::ResizeTo( int height )
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
    values_.resize( localHeight );
}

template<typename F>
const DistVector<F>& 
DistVector<F>::operator=( const DistVector<F>& x )
{
#ifndef RELEASE
    PushCallStack("DistVector::operator=");
#endif
    height_ = x.height_;
    SetComm( x.comm_ );
    values_ = x.values_;
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

} // namespace cliq

#endif // CLIQUE_DIST_VECTOR_MAIN_HPP
