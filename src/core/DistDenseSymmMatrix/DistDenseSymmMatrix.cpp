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

#include "./Chol-incl.hpp"
#include "./LDL-incl.hpp"

template<typename F>
clique::DistDenseSymmMatrix<F>::DistDenseSymmMatrix
( mpi::Comm comm, int gridHeight, int gridWidth )
: height_(0), blockSize_(1), 
  comm_(comm), gridHeight_(gridHeight), gridWidth_(gridWidth)
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::DistDenseSymmMatrix");
#endif
    const int commSize = mpi::CommSize( comm );
    if( commSize != gridHeight*gridWidth )
        throw std::logic_error
        ("Comm size was not compatible with specified grid dimensions");

    // Create a cartesian communicator
    int dimensions[2] = { gridWidth, gridHeight };
    int periods[2] = { true, true };
    int reorder = false;
    mpi::CartCreate
    ( comm_, 2, dimensions, periods, reorder, cartComm_ );

    // Set up the column and row communicators
    int remainingDimensions[2];
    remainingDimensions[0] = false;
    remainingDimensions[1] = true;
    mpi::CartSub( cartComm_, remainingDimensions, colComm_ );
    remainingDimensions[0] = true;
    remainingDimensions[1] = false;
    mpi::CartSub( cartComm_, remainingDimensions, rowComm_ );
    gridRow_ = mpi::CommRank( colComm_ );
    gridCol_ = mpi::CommRank( rowComm_ );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
clique::DistDenseSymmMatrix<F>::DistDenseSymmMatrix
( int height, int blockSize, mpi::Comm comm, int gridHeight, int gridWidth )
: height_(height), blockSize_(blockSize), 
  comm_(comm), gridHeight_(gridHeight), gridWidth_(gridWidth)
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::DistDenseSymmMatrix");
#endif
    const int commSize = mpi::CommSize( comm );
    if( commSize != gridHeight*gridWidth )
        throw std::logic_error
        ("Comm size was not compatible with specified grid dimensions");
    
    // Create a cartesian communicator
    int dimensions[2] = { gridWidth, gridHeight };
    int periods[2] = { true, true };
    int reorder = false;
    mpi::CartCreate
    ( comm_, 2, dimensions, periods, reorder, cartComm_ );

    // Set up the column and row communicators
    int remainingDimensions[2];
    remainingDimensions[0] = false;
    remainingDimensions[1] = true;
    mpi::CartSub( cartComm_, remainingDimensions, colComm_ );
    remainingDimensions[0] = true;
    remainingDimensions[1] = false;
    mpi::CartSub( cartComm_, remainingDimensions, rowComm_ );
    gridRow_ = mpi::CommRank( colComm_ );
    gridCol_ = mpi::CommRank( rowComm_ );

    // Determine the local height information
    const int localHeight = 
        LocalLength( height, blockSize, 0, gridRow_, gridHeight_ );
    const int localBlockHeight = BlockLength( localHeight, blockSize );

    // Determine the local width information
    const int localWidth = 
        LocalLength( height, blockSize, 0, gridCol_, gridWidth_ );
    const int localBlockWidth = BlockLength( localWidth, blockSize );

    // Compute the offsets for block column-major packed storage
    int totalLocalSize = 0;
    blockColHeights_.resize( localBlockWidth );
    blockColWidths_.resize( localBlockWidth );
    blockColRowOffsets_.resize( localBlockWidth );
    blockColColOffsets_.resize( localBlockWidth );
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        const int jBlock = gridCol_ + jLocalBlock*gridWidth;
        const int iLocalBlock = (jBlock-gridRow_+gridHeight-1)/gridHeight;
        const int iBlock = gridRow_ + iLocalBlock*gridHeight;
        blockColRowOffsets_[jLocalBlock] = iBlock*blockSize_;
        blockColColOffsets_[jLocalBlock] = jBlock*blockSize_;

        const int remainingLocalHeight = localHeight - iLocalBlock*blockSize;
        const int remainingLocalWidth = localWidth - jLocalBlock*blockSize;
        const int thisLocalHeight = remainingLocalHeight;
        const int thisLocalWidth = std::min(remainingLocalWidth,blockSize);
        blockColHeights_[jLocalBlock] = thisLocalHeight;
        blockColWidths_[jLocalBlock] = thisLocalWidth;
        totalLocalSize += thisLocalHeight*thisLocalWidth;
    }

    // Create space for the storage and fill in pointers to the block columns
    buffer_.resize( totalLocalSize );
    blockColBuffers_.resize( localBlockWidth );
    totalLocalSize = 0;
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        blockColBuffers_[jLocalBlock] = &buffer_[totalLocalSize];
        totalLocalSize += 
            blockColHeights_[jLocalBlock]*blockColWidths_[jLocalBlock];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::Reconfigure( int height, int blockSize )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::Reconfigure");
#endif
    height_ = height;
    blockSize_ = blockSize;

    // Clear the old storage
    buffer_.clear();
    blockColBuffers_.clear();
    blockColHeights_.clear();
    blockColWidths_.clear();

    // Determine the local height information
    const int localHeight = 
        LocalLength( height, blockSize, 0, gridRow_, gridHeight_ );
    const int localBlockHeight = BlockLength( localHeight, blockSize );

    // Determine the local width information
    const int localWidth = 
        LocalLength( height, blockSize, 0, gridCol_, gridWidth_ );
    const int localBlockWidth = BlockLength( localWidth, blockSize );

    // Compute the offsets for block column-major packed storage
    int totalLocalSize = 0;
    blockColHeights_.resize( localBlockWidth );
    blockColWidths_.resize( localBlockWidth );
    blockColRowOffsets_.resize( localBlockWidth );
    blockColColOffsets_.resize( localBlockWidth );
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        const int jBlock = gridCol_ + jLocalBlock*gridWidth_;
        const int iLocalBlock = (jBlock-gridRow_+gridHeight_-1)/gridHeight_;
        const int iBlock = gridRow_ + iLocalBlock*gridHeight_;
        blockColRowOffsets_[jLocalBlock] = iBlock*blockSize_;
        blockColColOffsets_[jLocalBlock] = jBlock*blockSize_;

        const int remainingLocalHeight = localHeight - iLocalBlock*blockSize;
        const int remainingLocalWidth = localWidth - jLocalBlock*blockSize;
        const int thisLocalHeight = remainingLocalHeight;
        const int thisLocalWidth = std::min(remainingLocalWidth,blockSize);
        blockColHeights_[jLocalBlock] = thisLocalHeight;
        blockColWidths_[jLocalBlock] = thisLocalWidth;
        totalLocalSize += thisLocalHeight*thisLocalWidth;
    }

    // Create space for the storage and fill in pointers to the block columns
    buffer_.resize( totalLocalSize );
    blockColBuffers_.resize( localBlockWidth );
    totalLocalSize = 0;
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        blockColBuffers_[jLocalBlock] = &buffer_[totalLocalSize];
        totalLocalSize += 
            blockColHeights_[jLocalBlock]*blockColWidths_[jLocalBlock];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::Print( std::string s ) const
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::Print");
#endif
    const int commSize = mpi::CommSize( comm_ );
    const int commRank = mpi::CommRank( comm_ );

    std::vector<F> sendBuf( height_*height_, 0 );

    // Pack our local matrix
    const int numLocalBlockCols = blockColHeights_.size();
    for( int jLocalBlock=0; jLocalBlock<numLocalBlockCols; ++jLocalBlock )
    {
        const F* blockCol = blockColBuffers_[jLocalBlock];
        const int iOffset = blockColRowOffsets_[jLocalBlock];
        const int jOffset = blockColColOffsets_[jLocalBlock];
        const int blockColHeight = blockColHeights_[jLocalBlock];
        const int blockColWidth = blockColWidths_[jLocalBlock];
        const int numLocalBlocks = (blockColHeight+blockSize_-1)/blockSize_;

        const int lastBlockHeight = LastBlockSize( blockColHeight, blockSize_ );

        for( int y=0; y<blockColWidth; ++y )
        {
            int block = 0;
            int blockOffset = iOffset;
            for( ; block<numLocalBlocks-1; ++block )
            {
                const F* col = &blockCol[block*blockSize_+y*blockColHeight];
                F* sendCol = &sendBuf[blockOffset+(jOffset+y)*height_];
                std::memcpy( sendCol, col, blockSize_*sizeof(F) );
                blockOffset += blockSize_;
            }

            const F* col = &blockCol[block*blockSize_+y*blockColHeight];
            F* sendCol = &sendBuf[blockOffset+(jOffset+y)*height_];
            std::memcpy( sendCol, col, lastBlockHeight*sizeof(F) );
        }
    }

    // if we are the root, allocate a receive buffer
    std::vector<F> recvBuf;
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

template<typename F>
void
clique::DistDenseSymmMatrix<F>::MakeZero()
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::MakeZero");
#endif
    std::memset( &buffer_[0], 0, buffer_.size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::MakeIdentity()
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::MakeIdentity");
#endif
    std::memset( &buffer_[0], 0, buffer_.size() );
    const int localBlockWidth = blockColBuffers_.size();
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        // Check if we own the diagonal block, if so, set to identity
        const int rowOffset = blockColRowOffsets_[jLocalBlock];
        const int colOffset = blockColColOffsets_[jLocalBlock];
        if( rowOffset == colOffset )
        {
            const int blockColHeight = blockColHeights_[jLocalBlock];
            const int blockColWidth = blockColWidths_[jLocalBlock];
            F* blockCol = blockColBuffers_[jLocalBlock];
            for( int j=0; j<std::min(blockColHeight,blockColWidth); ++j )
                blockCol[j+j*blockColHeight] = (F)1; 
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template class clique::DistDenseSymmMatrix<float>;
template class clique::DistDenseSymmMatrix<double>;
template class clique::DistDenseSymmMatrix<std::complex<float> >;
template class clique::DistDenseSymmMatrix<std::complex<double> >;
