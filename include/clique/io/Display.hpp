/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

void Display( const Graph& graph, std::string title="Graph" );
void Display( const DistGraph& graph, std::string title="DistGraph" );

template<typename T>
void Display( const SparseMatrix<T>& A, std::string title="SparseMatrix" );
template<typename T>
void Display
( const SparseMatrix<Complex<T> >& A, std::string title="SparseMatrix" );

template<typename T>
void Display
( const DistSparseMatrix<T>& A, 
  std::string title="DistSparseMatrix" );
template<typename T>
void Display
( const DistSparseMatrix<Complex<T> >& A, 
  std::string title="DistSparseMatrix" );

//
// Implementation begins here
//

inline void
Display( const Graph& graph, std::string title )
{
#ifndef RELEASE
    CallStackEntry cse("Display [Graph]");
#endif
#ifdef HAVE_QT5
    Matrix<int>* graphMat = new Matrix<int>;
    const int m = graph.NumTargets();
    const int n = graph.NumSources();
    elem::Zeros( *graphMat, m, n );

    const int numEdges = graph.NumEdges();
    const int* srcBuf = graph.LockedSourceBuffer();
    const int* tgtBuf = graph.LockedTargetBuffer();
    for( int e=0; e<numEdges; ++e )
        graphMat->Set( tgtBuf[e], srcBuf[e], 1 );

    QString qTitle = QString::fromStdString( title );
    elem::SpyWindow* spyWindow = new elem::SpyWindow;
    spyWindow->Spy( graphMat, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    elem::ProcessEvents( 200 );
#else
    Print( graph, title );
#endif
}

inline void
Display( const DistGraph& graph, std::string title ) 
{
#ifndef RELEASE
    CallStackEntry cse("Display [DistGraph]");
#endif
#ifdef HAVE_QT5
    const mpi::Comm comm = graph.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    const int numLocalEdges = graph.NumLocalEdges();
    std::vector<int> edgeSizes(commSize), edgeOffsets(commSize);
    mpi::AllGather( &numLocalEdges, 1, &edgeSizes[0], 1, comm );
    int numEdges=0;
    for( int q=0; q<commSize; ++q )
    {
        edgeOffsets[q] = numEdges;
        numEdges += edgeSizes[q];
    }

    std::vector<int> sources, targets;
    if( commRank == 0 )
    {
        sources.resize( numEdges );
        targets.resize( numEdges );
    }
    mpi::Gather
    ( graph.LockedSourceBuffer(), numLocalEdges,
      &sources[0], &edgeSizes[0], &edgeOffsets[0], 0, comm );
    mpi::Gather
    ( graph.LockedTargetBuffer(), numLocalEdges,
      &targets[0], &edgeSizes[0], &edgeOffsets[0], 0, comm );

    if( commRank == 0 )
    {
        Matrix<int>* graphMat = new Matrix<int>;
        const int m = graph.NumTargets();
        const int n = graph.NumSources();
        elem::Zeros( *graphMat, m, n );
        for( int e=0; e<numEdges; ++e )
            graphMat->Set( targets[e], sources[e], 1 );

        QString qTitle = QString::fromStdString( title );
        elem::SpyWindow* spyWindow = new elem::SpyWindow;
        spyWindow->Spy( graphMat, qTitle );
        spyWindow->show();

        // Spend at most 200 milliseconds rendering
        elem::ProcessEvents( 200 );
    }
#else
    Print( graph, title );
#endif
}

template<typename T>
inline void
Display( const SparseMatrix<T>& A, std::string title )
{
#ifndef RELEASE
    CallStackEntry cse("Print [SparseMatrix]");
#endif
#ifdef HAVE_QT5
    Matrix<double>* AFull = new Matrix<double>;
    const int m = A.Height();
    const int n = A.Width();
    elem::Zeros( *AFull, m, n );

    const int numEntries = A.NumEntries();
    const int* srcBuf = A.LockedSourceBuffer();
    const int* tgtBuf = A.LockedTargetBuffer();
    const T* valBuf = A.LockedValueBuffer();
    for( int s=0; s<numEntries; ++s )
        AFull->Set( tgtBuf[s], srcBuf[s], double(valBuf[s]) );

    QString qTitle = QString::fromStdString( title );
    elem::DisplayWindow* displayWindow = new elem::DisplayWindow;
    displayWindow->Display( AFull, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    elem::ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename T>
inline void
Display( const SparseMatrix<Complex<T> >& A, std::string title )
{
#ifndef RELEASE
    CallStackEntry cse("Print [SparseMatrix]");
#endif
#ifdef HAVE_QT5
    Matrix<Complex<double> >* AFull = new Matrix<Complex<double> >;
    const int m = A.Height();
    const int n = A.Width();
    elem::Zeros( *AFull, m, n );

    const int numEntries = A.NumEntries();
    const int* srcBuf = A.LockedSourceBuffer();
    const int* tgtBuf = A.LockedTargetBuffer();
    const Complex<T>* valBuf = A.LockedValueBuffer();
    for( int s=0; s<numEntries; ++s )
    {
        const Complex<double> alpha = 
            Complex<double>(valBuf[s].real,valBuf[s].imag);
        AFull->Set( tgtBuf[s], srcBuf[s], alpha );
    }

    QString qTitle = QString::fromStdString( title );
    elem::ComplexDisplayWindow* displayWindow = new elem::ComplexDisplayWindow;
    displayWindow->Display( AFull, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    elem::ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename T>
inline void
Display( const DistSparseMatrix<T>& A, std::string title )
{
#ifndef RELEASE
    CallStackEntry cse("Display [DistSparseMatrix]");
#endif
#ifdef HAVE_QT5
    const mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    const int numLocalEntries = A.NumLocalEntries();
    std::vector<int> nonzeroSizes(commSize), nonzeroOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &nonzeroSizes[0], 1, comm );
    int numNonzeros=0;
    for( int q=0; q<commSize; ++q )
    {
        nonzeroOffsets[q] = numNonzeros;
        numNonzeros += nonzeroSizes[q];
    }

    std::vector<int> sources, targets;
    std::vector<T> values;
    if( commRank == 0 )
    {
        sources.resize( numNonzeros );
        targets.resize( numNonzeros );
        values.resize( numNonzeros );
    }
    mpi::Gather
    ( A.LockedSourceBuffer(), numLocalEntries,
      &sources[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );
    mpi::Gather
    ( A.LockedTargetBuffer(), numLocalEntries,
      &targets[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );
    mpi::Gather
    ( A.LockedValueBuffer(), numLocalEntries,
      &values[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );

    if( commRank == 0 )
    {
        Matrix<double>* AFull = new Matrix<double>;
        const int m = A.Height();
        const int n = A.Width();
        elem::Zeros( *AFull, m, n );

        for( int s=0; s<numNonzeros; ++s )
            AFull->Set( targets[s], sources[s], double(values[s]) );

        QString qTitle = QString::fromStdString( title );
        elem::DisplayWindow* displayWindow = new elem::DisplayWindow;
        displayWindow->Display( AFull, qTitle );
        displayWindow->show();

        // Spend at most 200 milliseconds rendering
        elem::ProcessEvents( 200 );
    }
#else
    Print( A, title );
#endif
}

template<typename T>
inline void
Display( const DistSparseMatrix<Complex<T> >& A, std::string title )
{
#ifndef RELEASE
    CallStackEntry cse("Display [DistSparseMatrix]");
#endif
#ifdef HAVE_QT5
    const mpi::Comm comm = A.Comm();
    const int commSize = mpi::CommSize( comm );
    const int commRank = mpi::CommRank( comm );

    const int numLocalEntries = A.NumLocalEntries();
    std::vector<int> nonzeroSizes(commSize), nonzeroOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &nonzeroSizes[0], 1, comm );
    int numNonzeros=0;
    for( int q=0; q<commSize; ++q )
    {
        nonzeroOffsets[q] = numNonzeros;
        numNonzeros += nonzeroSizes[q];
    }

    std::vector<int> sources, targets;
    std::vector<Complex<T> > values;
    if( commRank == 0 )
    {
        sources.resize( numNonzeros );
        targets.resize( numNonzeros );
        values.resize( numNonzeros );
    }
    mpi::Gather
    ( A.LockedSourceBuffer(), numLocalEntries,
      &sources[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );
    mpi::Gather
    ( A.LockedTargetBuffer(), numLocalEntries,
      &targets[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );
    mpi::Gather
    ( A.LockedValueBuffer(), numLocalEntries,
      &values[0], &nonzeroSizes[0], &nonzeroOffsets[0], 0, comm );

    if( commRank == 0 )
    {
        Matrix<Complex<double> >* AFull = new Matrix<Complex<double> >;
        const int m = A.Height();
        const int n = A.Width();
        elem::Zeros( *AFull, m, n );

        for( int s=0; s<numNonzeros; ++s )
        {
            const Complex<double> alpha = 
                Complex<double>(values[s].real,values[s].imag);
            AFull->Set( targets[s], sources[s], alpha );
        }

        QString qTitle = QString::fromStdString( title );
        elem::ComplexDisplayWindow* displayWindow = 
            new elem::ComplexDisplayWindow;
        displayWindow->Display( AFull, qTitle );
        displayWindow->show();

        // Spend at most 200 milliseconds rendering
        elem::ProcessEvents( 200 );
    }
#else
    Print( A, title );
#endif
}

} // namespace cliq
