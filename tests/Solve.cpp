/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "clique.hpp"
using namespace cliq;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n1 = Input("--n1","first grid dimension",30);
        const int n2 = Input("--n2","second grid dimension",30);
        const int n3 = Input("--n3","third grid dimension",30);
        const int numRhs = Input("--numRhs","number of right-hand sides",5);
        const bool solve2d = Input("--solve2d","use 2d solve?",false);
        const bool selInv = Input("--selInv","selectively invert?",false);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        const int N = n1*n2*n3;
        DistSparseMatrix<double> A( N, comm );

        // Fill our portion of the 3D negative Laplacian using a n1 x n2 x n3
        // 7-point stencil in natural ordering: (x,y,z) at x + y*n1 + z*n1*n2
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of matrix...";
            std::cout.flush();
        }
        const double fillStart = mpi::Time();
        const int firstLocalRow = A.FirstLocalRow();
        const int localHeight = A.LocalHeight();
        A.StartAssembly();
        A.Reserve( 7*localHeight );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = firstLocalRow + iLocal;
            const int x = i % n1;
            const int y = (i/n1) % n2;
            const int z = i/(n1*n2);

            A.Update( i, i, 6. );
            if( x != 0 )
                A.Update( i, i-1, -1. );
            if( x != n1-1 )
                A.Update( i, i+1, -1. );
            if( y != 0 )
                A.Update( i, i-n1, -1. );
            if( y != n2-1 )
                A.Update( i, i+n1, -1. );
            if( z != 0 )
                A.Update( i, i-n1*n2, -1. );
            if( z != n3-1 )
                A.Update( i, i+n1*n2, -1. );
        } 
        A.StopAssembly();
        mpi::Barrier( comm );
        const double fillStop =  mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << fillStop-fillStart << " seconds" 
                      << std::endl;
        if( display )
        {
            Display( A );
            Display( A.DistGraph() );
        }
        if( print )
        {
            Print( A );
            Print( A.DistGraph() );
        }

        if( commRank == 0 )
        {
            std::cout << "Generating random X and forming Y := A X...";
            std::cout.flush();
        }
        const double multiplyStart = mpi::Time();
        DistMultiVec<double> X( N, numRhs, comm ), Y( N, numRhs, comm );
        MakeUniform( X );
        MakeZeros( Y );
        Multiply( 1., A, X, 0., Y );
        std::vector<double> YOrigNorms;
        Norms( Y, YOrigNorms );
        mpi::Barrier( comm );
        const double multiplyStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << multiplyStop-multiplyStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running nested dissection...";
            std::cout.flush();
        }
        const double nestedStart = mpi::Time();
        const DistGraph& graph = A.DistGraph();
        DistSymmInfo info;
        DistSeparatorTree sepTree;
        DistMap map, inverseMap;
        NestedDissection
        ( graph, map, sepTree, info, 
          sequential, numDistSeps, numSeqSeps, cutoff );
        map.FormInverse( inverseMap );
        mpi::Barrier( comm );
        const double nestedStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << nestedStop-nestedStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            const int numDistNodes = info.distNodes.size();
            const int numLocalNodes = info.localNodes.size();
            const int rootSepSize = info.distNodes.back().size;
            std::cout << "\n"
                      << "On the root process:\n"
                      << "-----------------------------------------\n"
                      << numLocalNodes << " local nodes\n"
                      << numDistNodes  << " distributed nodes\n"
                      << rootSepSize << " vertices in root separator\n"
                      << std::endl;
        }
        if( display )
        {
            std::ostringstream osBefore, osAfter;
            osBefore << "Structure before fact. on process " << commRank;
            osAfter << "Structure after fact. on process " << commRank;
            DisplayLocal( info, false, osBefore.str() );
            DisplayLocal( info, true, osAfter.str() );
        }

        if( commRank == 0 )
        {
            std::cout << "Building DistSymmFrontTree...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double buildStart = mpi::Time();
        DistSymmFrontTree<double> frontTree( TRANSPOSE, A, map, sepTree, info );
        mpi::Barrier( comm );
        const double buildStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << buildStop-buildStart << " seconds"
                      << std::endl;

        double numLocalEntries, minLocalEntries, maxLocalEntries, 
               numGlobalEntries;
        frontTree.MemoryInfo
        ( numLocalEntries, minLocalEntries, maxLocalEntries, numGlobalEntries );
        double numLocalFlops, minLocalFlops, maxLocalFlops, numGlobalFlops;
        frontTree.FactorizationWork
        ( numLocalFlops, minLocalFlops, maxLocalFlops, numGlobalFlops );
        if( commRank == 0 )
        {
            std::cout << "Original memory usage for fronts...\n"
                      << "  min local: " << minLocalEntries*sizeof(double)/1e6 
                      << " MB\n"
                      << "  max local: " << maxLocalEntries*sizeof(double)/1e6 
                      << " MB\n"
                      << "  global:    " << numGlobalEntries*sizeof(double)/1e6
                      << " MB\n"
                      << "\n"
                      << "Factorization work...\n"
                      << "  min local: " << minLocalFlops/1.e9 << " GFlops\n"
                      << "  max local: " << maxLocalFlops/1.e9 << " GFlops\n"
                      << "  global:    " << numGlobalFlops/1.e9 << " GFlops\n"
                      << std::endl;
        }

        if( commRank == 0 )
        {
            std::cout << "Running LDL^T and redistribution...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double ldlStart = mpi::Time();
        if( solve2d )
        {
            if( selInv )
                LDL( info, frontTree, LDL_SELINV_2D );
            else
                LDL( info, frontTree, LDL_2D );
        }
        else
        {
            if( selInv )
                LDL( info, frontTree, LDL_SELINV_1D );
            else
                LDL( info, frontTree, LDL_1D );
        }
        mpi::Barrier( comm );
        const double ldlStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << ldlStop-ldlStart << " seconds" 
                      << std::endl;

        if( commRank == 0 )
            std::cout << "Memory usage for fronts after factorization..."
                      << std::endl;
        frontTree.MemoryInfo
        ( numLocalEntries, minLocalEntries, maxLocalEntries, numGlobalEntries );
        if( commRank == 0 )
        {
            std::cout << "  min local: " << minLocalEntries*sizeof(double)/1e6 
                      << " MB\n"
                      << "  max local: " << maxLocalEntries*sizeof(double)/1e6 
                      << " MB\n"
                      << "  global:    " << numGlobalEntries*sizeof(double)/1e6
                      << " MB\n"
                      << std::endl;
        }

        if( commRank == 0 )
        {
            std::cout << "Solving against Y...";
            std::cout.flush();
        }
        const double solveStart = mpi::Time();
        if( solve2d )
        {
            DistNodalMatrix<double> YNodal;
            YNodal.Pull( inverseMap, info, Y );
            Solve( info, frontTree, YNodal );
            YNodal.Push( inverseMap, info, Y );
        }
        else
        {
            DistNodalMultiVec<double> YNodal;
            YNodal.Pull( inverseMap, info, Y );
            Solve( info, frontTree, YNodal );
            YNodal.Push( inverseMap, info, Y );
        }
        mpi::Barrier( comm );
        const double solveStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << solveStop-solveStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
            std::cout << "Checking error in computed solution..." << std::endl;
        std::vector<double> XNorms, YNorms;
        Norms( X, XNorms );
        Norms( Y, YNorms );
        Axpy( -1., X, Y );
        std::vector<double> errorNorms;
        Norms( Y, errorNorms );
        if( commRank == 0 )
        {
            for( int j=0; j<numRhs; ++j )
            {
                std::cout << "Right-hand side " << j << ":\n"
                          << "|| x     ||_2 = " << XNorms[j] << "\n"
                          << "|| xComp ||_2 = " << YNorms[j] << "\n"
                          << "|| A x   ||_2 = " << YOrigNorms[j] << "\n"
                          << "|| error ||_2 = " << errorNorms[j] << "\n"
                          << std::endl;
            }
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
