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

// y := alpha A x + beta y
template<typename T>
void Multiply
( T alpha, const DistSparseMatrix<T>& A, const DistVector<T>& x,
  T beta,                                      DistVector<T>& y );

// Y := alpha A X + beta Y
template<typename T>
void Multiply
( T alpha, const DistSparseMatrix<T>& A, const DistMultiVector<T>& X,
  T beta,                                      DistMultiVector<T>& Y );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
void Multiply
( T alpha, const DistSparseMatrix<T>& A, const DistVector<T>& x,
  T beta,                                      DistVector<T>& y )
{
#ifndef RELEASE
    PushCallStack("Multiply");
    if( A.Height() != y.Height() || A.Width() != x.Height() )
        throw std::logic_error("A, x, and y did not conform");
    if( !mpi::CongruentComms( A.Comm(), x.Comm() ) || 
        !mpi::CongruentComms( x.Comm(), y.Comm() ) )
        throw std::logic_error("Communicators did not match");
#endif
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    const int xLocalHeight = x.LocalHeight();
    const int yLocalHeight = y.LocalHeight();
    const int blocksize = A.Blocksize();
    const int firstLocalRow = A.FirstLocalRow();
    const int numLocalEntries = A.NumLocalEntries();

    // y := beta y
    for( int iLocal=0; iLocal<yLocalHeight; ++iLocal )
        y.SetLocal( iLocal, beta*y.GetLocal(iLocal) );

    // Compute the set of row indices that we need from x
    std::set<int> indexSet;
    for( int e=0; e<numLocalEntries; ++e )
        indexSet.insert( A.Col(e) );
    const int numRecvIndices = indexSet.size();
    std::vector<int> recvIndices( numRecvIndices );
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
            recvIndices[offset++] = j;
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
    int numSendIndices=0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSendIndices;
        numSendIndices += sendSizes[q];
    }
    std::vector<int> sendIndices( numSendIndices );
    mpi::AllToAll
    ( &recvIndices[0], &recvSizes[0], &recvOffsets[0],
      &sendIndices[0], &sendSizes[0], &sendOffsets[0], comm );

    // Pack the send values
    std::vector<T> sendValues( numSendIndices );
    for( int s=0; s<numSendIndices; ++s )
    {
        const int i = sendIndices[s];
        const int iLocal = i - firstLocalRow;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= xLocalHeight )
            throw std::logic_error("iLocal was out of bounds");
#endif
        sendValues[s] = x.GetLocal( iLocal );
    }

    // Send them back
    std::vector<T> recvValues( numRecvIndices );
    mpi::AllToAll
    ( &sendValues[0], &sendSizes[0], &sendOffsets[0],
      &recvValues[0], &recvSizes[0], &recvOffsets[0], comm );
     
    // Perform the local multiply-accumulate, y := alpha A x + y
    std::vector<int>::const_iterator vecIt;
    for( int iLocal=0; iLocal<yLocalHeight; ++iLocal )
    {
        const int offset = A.LocalEntryOffset( iLocal );
        const int rowSize = A.NumConnections( iLocal );
        for( int k=0; k<rowSize; ++k )
        {
            const int col = A.Col( offset+k );
            vecIt = std::lower_bound
            ( recvIndices.begin(), recvIndices.end(), col );
#ifndef RELEASE            
            if( vecIt == recvIndices.end() )
                throw std::logic_error("Could not find recv index");
#endif
            const int colOffset = vecIt - recvIndices.begin();
            const T AValue = A.Value(k+offset);
            const T xValue = recvValues[colOffset];
            const T update = alpha*AValue*xValue;
            y.UpdateLocal( iLocal, update );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void Multiply
( T alpha, const DistSparseMatrix<T>& A, const DistMultiVector<T>& X,
  T beta,                                      DistMultiVector<T>& Y )
{
#ifndef RELEASE
    PushCallStack("Multiply");
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
    const int numRecvIndices = indexSet.size();
    std::vector<int> recvIndices( numRecvIndices );
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
            recvIndices[offset++] = j;
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
    int numSendIndices=0;
    std::vector<int> sendOffsets( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = numSendIndices;
        numSendIndices += sendSizes[q];
    }
    std::vector<int> sendIndices( numSendIndices );
    mpi::AllToAll
    ( &recvIndices[0], &recvSizes[0], &recvOffsets[0],
      &sendIndices[0], &sendSizes[0], &sendOffsets[0], comm );

    // Pack the send values
    std::vector<T> sendValues( numSendIndices*width );
    for( int s=0; s<numSendIndices; ++s )
    {
        const int i = sendIndices[s];
        const int iLocal = i - firstLocalRow;
#ifndef RELEASE
        if( iLocal < 0 || iLocal >= XLocalHeight )
            throw std::logic_error("iLocal was out of bounds");
#endif
        for( int j=0; j<width; ++j )
            sendValues[s*width+j] = X.GetLocal( iLocal, j );
    }

    // Send them back
    std::vector<T> recvValues( numRecvIndices*width );
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
    std::vector<int>::const_iterator vecIt;
    for( int iLocal=0; iLocal<YLocalHeight; ++iLocal )
    {
        const int offset = A.LocalEntryOffset( iLocal );
        const int rowSize = A.NumConnections( iLocal );
        for( int k=0; k<rowSize; ++k )
        {
            const int col = A.Col( offset+k );
            vecIt = std::lower_bound
            ( recvIndices.begin(), recvIndices.end(), col );
#ifndef RELEASE            
            if( vecIt == recvIndices.end() )
                throw std::logic_error("Could not find recv index");
#endif
            const int colOffset = vecIt - recvIndices.begin();
            const T AValue = A.Value(k+offset);
            for( int j=0; j<width; ++j )
            {
                const T XValue = recvValues[colOffset*width+j];
                const T update = alpha*AValue*XValue;
                Y.UpdateLocal( iLocal, j, update );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
