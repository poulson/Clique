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
( int height, int blockSize, MPI_Comm comm, int gridHeight, int gridWidth )
: height_(height), blockSize_(blockSize), 
  comm_(comm), gridHeight_(gridHeight), gridWidth_(gridWidth)
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::DistDenseSymmMatrix");
#endif
    int commSize, commRank;
    MPI_Comm_size( comm, &commSize );
    MPI_Comm_rank( comm, &commRank );

    if( commSize != gridHeight*gridWidth )
        throw std::logic_error
        ("Comm size was not compatible with specified grid dimensions");
    gridRow_ = commRank % gridHeight_;
    gridCol_ = commRank / gridHeight_;

    // Determine the local height information
    const int offsetHeight = std::max(0,height-gridRow_*blockSize);
    const int remainderHeight = offsetHeight % (gridHeight*blockSize);
    const int localHeight = (offsetHeight/(gridHeight*blockSize))*blockSize +
        std::min(blockSize,remainderHeight);
    const int localBlockHeight = 
        ( remainderHeight==0 ? 
          localHeight/blockSize : localHeight/blockSize+1 );

    // Determine the local width information
    const int offsetWidth = std::max(0,height-gridCol_*blockSize);
    const int remainderWidth = offsetWidth % (gridWidth*blockSize);
    const int localWidth = (offsetWidth/(gridWidth*blockSize))*blockSize +
        std::min(blockSize,remainderWidth);
    const int localBlockWidth =
        ( remainderWidth==0 ?
          localWidth/blockSize : localWidth/blockSize+1 );

    // Compute the offsets for block column-major packed storage
    int totalLocalSize = 0;
    blockColumnHeights_.resize( localBlockWidth );
    blockColumnWidths_.resize( localBlockWidth );
    for( int jLocalBlock=0; jLocalBlock<localBlockWidth; ++jLocalBlock )
    {
        const int jBlock = gridCol_ + jLocalBlock*gridWidth;
        const int iLocalBlock = (jBlock-gridRow_+gridHeight-1)/gridHeight;
        const int thisLocalHeight = localHeight-iLocalBlock*blockSize;
        const int thisLocalWidth = localWidth-jLocalBlock*blockSize;

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
