/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

// Y := alpha A X + beta Y
template<typename T>
void Multiply
( T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X,
  T beta,                                      DistMultiVec<T>& Y );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
void Multiply
( T alpha, const DistSparseMatrix<T>& A, const DistMultiVec<T>& X,
  T beta,                                      DistMultiVec<T>& Y )
{
#ifndef RELEASE
    CallStackEntry entry("Multiply");
    if( A.Height() != Y.Height() || A.Width() != X.Height() || 
        X.Width() != Y.Width() )
        throw std::logic_error("A, X, and Y did not conform");
    if( !mpi::CongruentComms( A.Comm(), X.Comm() ) || 
        !mpi::CongruentComms( X.Comm(), Y.Comm() ) )
        throw std::logic_error("Communicators did not match");
#endif
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    const int XLocalHeight = X.LocalHeight();
    const int YLocalHeight = Y.LocalHeight();
    const int width = X.Width();
    const int blocksize = A.Blocksize();
    const int firstLocalRow = A.FirstLocalRow();
    const int numLocalEntries = A.NumLocalEntries();

    // Y := beta Y
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<YLocalHeight; ++iLocal )
            Y.SetLocal( iLocal, j, beta*Y.GetLocal(iLocal,j) );

    // Compute the set of row indices that we need from X
    std::set<int> indexSet;
    for( int e=0; e<numLocalEntries; ++e )
        indexSet.insert( A.Col(e) );
    const int numRecvInd = indexSet.size();
    std::vector<int> recvInd( numRecvInd );
    std::vector<int> recvSizes( commSize, 0 );
    std::vector<int> recvOffsets( commSize );
    {
        int offset=0, lastOffset=0, qPrev=0;
        std::set<int>::const_iterator setIt;
        for( setIt=indexSet.begin(); setIt!=indexSet.end(); ++setIt )
        {
            const int j = *setIt;
            const int q = RowToProcess( j, blocksize, commSize );
            while( qPrev != q )
            {
                recvSizes[qPrev] = offset - lastOffset;
                recvOffsets[qPrev+1] = offset;

                lastOffset = offset;
                ++qPrev;
            }
            recvInd[offset++] = j;
        }
        while( qPrev != commSize-1 )
        {
            recvSizes[qPrev] = offset - lastOffset;
            recvOffsets[qPrev+1] = offset;
            lastOffset = offset;
            ++qPrev;
        }
        recvSizes[commSize-1] = offset - lastOffset;
    }

    // Coordinate
    std::vector<int> sendSizes( commSize );
    mpi::AllToAll( &recvSizes[0], 1, &sendSizes[0], 1, comm );
    int numSendInd=0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSendInd;
        numSendInd += sendSizes[q];
    }
    std::vector<int> sendInd( numSendInd );
    mpi::AllToAll
    ( &recvInd[0], &recvSizes[0], &recvOffsets[0],
      &sendInd[0], &sendSizes[0], &sendOffsets[0], comm );

    // Pack the send values
    std::vector<T> sendValues( numSendInd*width );
    for( int s=0; s<numSendInd; ++s )
    {
        const int i = sendInd[s];
        const int iLocal = i - firstLocalRow;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= XLocalHeight )
            throw std::logic_error("iLocal was out of bounds");
#endif
        for( int j=0; j<width; ++j )
            sendValues[s*width+j] = X.GetLocal( iLocal, j );
    }

    // Send them back
    std::vector<T> recvValues( numRecvInd*width );
    for( int q=0; q<commSize; ++q )
    {
        sendSizes[q] *= width;
        sendOffsets[q] *= width;
        recvSizes[q] *= width;
        recvOffsets[q] *= width;
    }
    mpi::AllToAll
    ( &sendValues[0], &sendSizes[0], &sendOffsets[0],
      &recvValues[0], &recvSizes[0], &recvOffsets[0], comm );
     
    // Perform the local multiply-accumulate, y := alpha A x + y
    for( int iLocal=0; iLocal<YLocalHeight; ++iLocal )
    {
        const int offset = A.LocalEntryOffset( iLocal );
        const int rowSize = A.NumConnections( iLocal );
        for( int k=0; k<rowSize; ++k )
        {
            const int col = A.Col( offset+k );
            const int colOffset = Find( recvInd, col );
            const T AValue = A.Value(k+offset);
            for( int j=0; j<width; ++j )
            {
                const T XValue = recvValues[colOffset*width+j];
                const T update = alpha*AValue*XValue;
                Y.UpdateLocal( iLocal, j, update );
            }
        }
    }
}

} // namespace cliq
