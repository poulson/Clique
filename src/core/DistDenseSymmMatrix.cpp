/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
#include "clique.hpp"

template<typename T>
clique::DistDenseSymmMatrix<T>::DistDenseSymmMatrix
( mpi::Comm comm, int gridHeight, int gridWidth )
: height_(0), blockSize_(1), 
  comm_(comm), gridHeight_(gridHeight), gridWidth_(gridWidth)
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::DistDenseSymmMatrix");
#endif
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    if( commSize != gridHeight*gridWidth )
        throw std::logic_error
        ("Comm size was not compatible with specified grid dimensions");
    gridRow_ = commRank % gridHeight_;
    gridCol_ = commRank / gridHeight_;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
clique::DistDenseSymmMatrix<T>::DistDenseSymmMatrix
( int height, int blockSize, mpi::Comm comm, int gridHeight, int gridWidth )
: height_(height), blockSize_(blockSize), 
  comm_(comm), gridHeight_(gridHeight), gridWidth_(gridWidth)
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::DistDenseSymmMatrix");
#endif
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    if( commSize != gridHeight*gridWidth )
        throw std::logic_error
        ("Comm size was not compatible with specified grid dimensions");
    gridRow_ = commRank % gridHeight_;
    gridCol_ = commRank / gridHeight_;

    // Determine the local height information
    const int offsetHeight = std::max(0,height-gridRow_*blockSize);
    const int remainderHeight = offsetHeight % blockSize;
    const int localHeight = (offsetHeight/(gridHeight*blockSize))*blockSize +
        std::min(blockSize,offsetHeight%(gridHeight*blockSize));
    const int localBlockHeight = 
        ( remainderHeight==0 ? 
          localHeight/blockSize : localHeight/blockSize+1 );

    // Determine the local width information
    const int offsetWidth = std::max(0,height-gridCol_*blockSize);
    const int remainderWidth = offsetWidth % blockSize;
    const int localWidth = (offsetWidth/(gridWidth*blockSize))*blockSize +
        std::min(blockSize,offsetWidth%(gridWidth*blockSize));
    const int localBlockWidth =
        ( remainderWidth==0 ?
          localWidth/blockSize : localWidth/blockSize+1 );

    // Compute the offsets for block column-major packed storage
    int totalLocalSize = 0;
    blockColumnHeights_.resize( localBlockWidth );
    blockColumnWidths_.resize( localBlockWidth );
    blockColumnRowOffsets_.resize( localBlockWidth );
    blockColumnColumnOffsets_.resize( localBlockWidth );
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        const int jBlock = gridCol_ + jLocalBlock*gridWidth;
        const int iLocalBlock = (jBlock-gridRow_+gridHeight-1)/gridHeight;
        const int iBlock = gridRow_ + iLocalBlock*gridHeight;
        blockColumnRowOffsets_[jLocalBlock] = iBlock*blockSize_;
        blockColumnColumnOffsets_[jLocalBlock] = jBlock*blockSize_;

        const int remainingLocalHeight = localHeight - iLocalBlock*blockSize;
        const int remainingLocalWidth = localWidth - jLocalBlock*blockSize;
        const int thisLocalHeight = remainingLocalHeight;
        const int thisLocalWidth = std::min(remainingLocalWidth,blockSize);
        blockColumnHeights_[jLocalBlock] = thisLocalHeight;
        blockColumnWidths_[jLocalBlock] = thisLocalWidth;
        totalLocalSize += thisLocalHeight*thisLocalWidth;
    }

    // Create space for the storage and fill in pointers to the block columns
    buffer_.resize( totalLocalSize );
    blockColumnBuffers_.resize( localBlockWidth );
    totalLocalSize = 0;
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        blockColumnBuffers_[jLocalBlock] = &buffer_[totalLocalSize];
        totalLocalSize += 
            blockColumnHeights_[jLocalBlock]*blockColumnWidths_[jLocalBlock];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
clique::DistDenseSymmMatrix<T>::Reconfigure( int height, int blockSize )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::Reconfigure");
#endif
    height_ = height;
    blockSize_ = blockSize;

    // Clear the old storage
    buffer_.clear();
    blockColumnBuffers_.clear();
    blockColumnHeights_.clear();
    blockColumnWidths_.clear();

    // Determine the local height information
    const int offsetHeight = std::max(0,height-gridRow_*blockSize);
    const int remainderHeight = offsetHeight % blockSize;
    const int localHeight = (offsetHeight/(gridHeight_*blockSize))*blockSize +
        std::min(blockSize,offsetHeight%(gridHeight_*blockSize));
    const int localBlockHeight = 
        ( remainderHeight==0 ? 
          localHeight/blockSize : localHeight/blockSize+1 );

    // Determine the local width information
    const int offsetWidth = std::max(0,height-gridCol_*blockSize);
    const int remainderWidth = offsetWidth % blockSize;
    const int localWidth = (offsetWidth/(gridWidth_*blockSize))*blockSize +
        std::min(blockSize,offsetWidth%(gridWidth_*blockSize));
    const int localBlockWidth =
        ( remainderWidth==0 ?
          localWidth/blockSize : localWidth/blockSize+1 );

    // Compute the offsets for block column-major packed storage
    int totalLocalSize = 0;
    blockColumnHeights_.resize( localBlockWidth );
    blockColumnWidths_.resize( localBlockWidth );
    blockColumnRowOffsets_.resize( localBlockWidth );
    blockColumnColumnOffsets_.resize( localBlockWidth );
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        const int jBlock = gridCol_ + jLocalBlock*gridWidth_;
        const int iLocalBlock = (jBlock-gridRow_+gridHeight_-1)/gridHeight_;
        const int iBlock = gridRow_ + iLocalBlock*gridHeight_;
        blockColumnRowOffsets_[jLocalBlock] = iBlock*blockSize_;
        blockColumnColumnOffsets_[jLocalBlock] = jBlock*blockSize_;

        const int remainingLocalHeight = localHeight - iLocalBlock*blockSize;
        const int remainingLocalWidth = localWidth - jLocalBlock*blockSize;
        const int thisLocalHeight = remainingLocalHeight;
        const int thisLocalWidth = std::min(remainingLocalWidth,blockSize);
        blockColumnHeights_[jLocalBlock] = thisLocalHeight;
        blockColumnWidths_[jLocalBlock] = thisLocalWidth;
        totalLocalSize += thisLocalHeight*thisLocalWidth;
    }

    // Create space for the storage and fill in pointers to the block columns
    buffer_.resize( totalLocalSize );
    blockColumnBuffers_.resize( localBlockWidth );
    totalLocalSize = 0;
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        blockColumnBuffers_[jLocalBlock] = &buffer_[totalLocalSize];
        totalLocalSize += 
            blockColumnHeights_[jLocalBlock]*blockColumnWidths_[jLocalBlock];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
clique::DistDenseSymmMatrix<T>::Print( std::string s ) const
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::Print");
#endif
    const int commSize = mpi::CommSize( comm_ );
    const int commRank = mpi::CommRank( comm_ );

    std::vector<T> sendBuf( height_*height_, 0 );

    // Pack our local matrix
    const int numLocalBlockColumns = blockColumnHeights_.size();
    for( int jLocalBlock=0; jLocalBlock<numLocalBlockColumns; ++jLocalBlock )
    {
        const T* blockColumn = blockColumnBuffers_[jLocalBlock];
        const int iOffset = blockColumnRowOffsets_[jLocalBlock];
        const int jOffset = blockColumnColumnOffsets_[jLocalBlock];
        const int blockColumnHeight = blockColumnHeights_[jLocalBlock];
        const int blockColumnWidth = blockColumnWidths_[jLocalBlock];
        const int numLocalBlocks = (blockColumnHeight+blockSize_-1)/blockSize_;

        const int remainder = blockColumnHeight % blockSize_;
        const int lastBlockHeight = ( remainder==0 ? blockSize_ : remainder );

        for( int y=0; y<blockColumnWidth; ++y )
        {
            int block = 0;
            int blockOffset = iOffset;
            for( ; block<numLocalBlocks-1; ++block )
            {
                const T* column = 
                    &blockColumn[block*blockSize_+y*blockColumnHeight];
                T* sendColumn = 
                    &sendBuf[blockOffset+(jOffset+y)*height_];
                std::memcpy( sendColumn, column, blockSize_*sizeof(T) );
                blockOffset += blockSize_;
            }

            if( numLocalBlocks != 0 )
            {
                const T* column = 
                    &blockColumn[block*blockSize_+y*blockColumnHeight];
                T* sendColumn = 
                    &sendBuf[blockOffset+(jOffset+y)*height_];
                std::memcpy( sendColumn, column, lastBlockHeight*sizeof(T) );
            }
        }
    }

    // if we are the root, allocate a receive buffer
    std::vector<T> recvBuf;
    if( commRank == 0 )
        recvBuf.resize( height_*height_ );

    // Sum the contributions and send to the root
    mpi::Reduce
    ( &sendBuf[0], &recvBuf[0], height_*height_, mpi::SUM, 0, comm_ );

    // Print the data from the root
    if( commRank == 0 )
    {
        if( s != "" )
            std::cout << s << "\n";
        for( int i=0; i<height_; ++i )
        {
            for( int j=0; j<height_; ++j )
                std::cout << recvBuf[i+j*height_] << " ";
            std::cout << "\n";
        }
        std::cout << std::endl;
    }
    mpi::Barrier( comm_ );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
clique::DistDenseSymmMatrix<T>::MakeZero()
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::MakeZero");
#endif
    std::memset( &buffer_[0], 0, buffer_.size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
clique::DistDenseSymmMatrix<T>::MakeIdentity()
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::MakeIdentity");
#endif
    std::memset( &buffer_[0], 0, buffer_.size() );
    const int localBlockWidth = blockColumnBuffers_.size();
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        // Check if we own the diagonal block, if so, set to identity
        const int rowOffset = blockColumnRowOffsets_[jLocalBlock];
        const int colOffset = blockColumnColumnOffsets_[jLocalBlock];
        if( rowOffset == colOffset )
        {
            const int blockColumnHeight = blockColumnHeights_[jLocalBlock];
            const int blockColumnWidth = blockColumnWidths_[jLocalBlock];
            T* blockColumn = blockColumnBuffers_[jLocalBlock];
            for( int j=0; j<std::min(blockColumnHeight,blockColumnWidth); ++j )
                blockColumn[j+j*blockColumnHeight] = (T)1; 
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template class clique::DistDenseSymmMatrix<int>;
template class clique::DistDenseSymmMatrix<float>;
template class clique::DistDenseSymmMatrix<double>;
template class clique::DistDenseSymmMatrix<std::complex<float> >;
template class clique::DistDenseSymmMatrix<std::complex<double> >;
