/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_IO_DISPLAY_HPP
#define CLIQ_IO_DISPLAY_HPP

namespace cliq {

void Display( const Graph& graph, std::string title="Graph" );
void Display( const DistGraph& graph, std::string title="DistGraph" );

template<typename T>
void Display( const SparseMatrix<T>& A, std::string title="SparseMatrix" );
template<typename T>
void Display
( const SparseMatrix<Complex<T>>& A, std::string title="SparseMatrix" );

template<typename T>
void Display
( const DistSparseMatrix<T>& A, 
  std::string title="DistSparseMatrix" );
template<typename T>
void Display
( const DistSparseMatrix<Complex<T>>& A, 
  std::string title="DistSparseMatrix" );

void DisplayLocal
( const DistSymmInfo& info, 
  bool beforeFact=true, std::string title="Local DistSymmInfo" );

//
// Implementation begins here
//

inline void
Display( const Graph& graph, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display [Graph]"))
#ifdef HAVE_QT5
    auto graphMat = new Matrix<int>;
    const int m = graph.NumTargets();
    const int n = graph.NumSources();
    elem::Zeros( *graphMat, m, n );

    const int numEdges = graph.NumEdges();
    const int* srcBuf = graph.LockedSourceBuffer();
    const int* tgtBuf = graph.LockedTargetBuffer();
    for( int e=0; e<numEdges; ++e )
        graphMat->Set( tgtBuf[e], srcBuf[e], 1 );

    QString qTitle = QString::fromStdString( title );
    auto spyWindow = new elem::SpyWindow;
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
    DEBUG_ONLY(CallStackEntry cse("Display [DistGraph]"))
#ifdef HAVE_QT5
    const mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

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
        auto graphMat = new Matrix<int>;
        const int m = graph.NumTargets();
        const int n = graph.NumSources();
        elem::Zeros( *graphMat, m, n );
        for( int e=0; e<numEdges; ++e )
            graphMat->Set( targets[e], sources[e], 1 );

        QString qTitle = QString::fromStdString( title );
        auto spyWindow = new elem::SpyWindow;
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
    DEBUG_ONLY(CallStackEntry cse("Print [SparseMatrix]"))
#ifdef HAVE_QT5
    auto AFull = new Matrix<double>;
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
    auto displayWindow = new elem::DisplayWindow;
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
Display( const SparseMatrix<Complex<T>>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Print [SparseMatrix]"))
#ifdef HAVE_QT5
    auto AFull = new Matrix<Complex<double>>;
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
    auto displayWindow = new elem::ComplexDisplayWindow;
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
    DEBUG_ONLY(CallStackEntry cse("Display [DistSparseMatrix]"))
#ifdef HAVE_QT5
    const mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

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
        auto AFull = new Matrix<double>;
        const int m = A.Height();
        const int n = A.Width();
        elem::Zeros( *AFull, m, n );

        for( int s=0; s<numNonzeros; ++s )
            AFull->Set( targets[s], sources[s], double(values[s]) );

        QString qTitle = QString::fromStdString( title );
        auto displayWindow = new elem::DisplayWindow;
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
Display( const DistSparseMatrix<Complex<T>>& A, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("Display [DistSparseMatrix]"))
#ifdef HAVE_QT5
    const mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

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
    std::vector<Complex<T>> values;
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
        auto AFull = new Matrix<Complex<double>>;
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
        auto displayWindow = new elem::ComplexDisplayWindow;
        displayWindow->Display( AFull, qTitle );
        displayWindow->show();

        // Spend at most 200 milliseconds rendering
        elem::ProcessEvents( 200 );
    }
#else
    Print( A, title );
#endif
}

inline void 
DisplayLocal( const DistSymmInfo& info, bool beforeFact, std::string title )
{
    DEBUG_ONLY(CallStackEntry cse("DisplayLocal [DistSymmInfo]"))
#ifdef HAVE_QT5
    const int n = info.distNodes.back().size + info.distNodes.back().off;
    auto graphMat = new Matrix<int>;
    elem::Zeros( *graphMat, n, n );

    const int numLocal = info.localNodes.size();
    for( int s=0; s<numLocal; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        for( int j=0; j<node.size; ++j )
            for( int i=0; i<node.size; ++i )
                graphMat->Set( i+node.off, j+node.off, 1 );
        if( beforeFact )
        {
            const int origStructSize = node.origLowerStruct.size();
            for( int i=0; i<origStructSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.origLowerStruct[i], j+node.off, 1 );
        }
        else
        {
            const int structSize = node.lowerStruct.size();
            for( int i=0; i<structSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.lowerStruct[i], j+node.off, 1 );
        }
    }

    const int numDist = info.distNodes.size();
    for( int s=0; s<numDist; ++s )
    {
        const DistSymmNodeInfo& node = info.distNodes[s];
        for( int j=0; j<node.size; ++j )
            for( int i=0; i<node.size; ++i )
                graphMat->Set( i+node.off, j+node.off, 1 );
        if( beforeFact )
        {
            const int origStructSize = node.origLowerStruct.size();
            for( int i=0; i<origStructSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.origLowerStruct[i], j+node.off, 1 );
        }
        else
        {
            const int structSize = node.lowerStruct.size();
            for( int i=0; i<structSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.lowerStruct[i], j+node.off, 1 );
        }
    }

    QString qTitle = QString::fromStdString( title );
    auto spyWindow = new elem::SpyWindow;
    spyWindow->Spy( graphMat, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    elem::ProcessEvents( 200 );
#else
    PrintLocal( info );
#endif
}

} // namespace cliq

#endif // ifndef CLIQ_IO_DISPLAY_HPP
