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
#include "clique.hpp"
using namespace cliq;

void Usage()
{
    std::cout
      << "DistSparseMatrix <n> [cutoff=128] [numDistSeps=10] [numSeqSeps=5]\n"
      << "  n: size of n x n x n mesh\n"
      << "  cutoff: maximum size of leaf node\n"
      << "  numDistSeps: number of distributed separators to try\n"
      << "  numSeqSeps: number of sequential separators to try\n"
      << std::endl;
}

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 2 )
    {
        if( commRank == 0 )
            Usage();
        cliq::Finalize();
        return 0;
    }
    const int n = atoi( argv[1] );
    const int cutoff = ( argc >= 3 ? atoi( argv[2] ) : 128 );
    const int numDistSeps = ( argc >= 4 ? atoi( argv[3] ) : 10 );
    const int numSeqSeps = ( argc >= 5 ? atoi( argv[4] ) : 5 );

    try
    {
        const int N = n*n*n;
        DistSparseMatrix<double> A( N, N, comm );

        const int firstLocalRow = A.FirstLocalRow();
        const int localHeight = A.LocalHeight();

        // Fill our portion of the 3D negative Laplacian using a n x n x n 
        // 7-point stencil in natural ordering: (x,y,z) at x + y*n + z*n*n
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of matrix...";
            std::cout.flush();
        }
        A.StartAssembly();
        A.Reserve( 7*localHeight );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = firstLocalRow + iLocal;
            const int x = i % n;
            const int y = (i/n) % n;
            const int z = i/(n*n);

            if( z != 0 )
                A.PushBack( i, i-n*n, -1. );
            if( y != 0 )
                A.PushBack( i, i-n, -1. );
            if( x != 0 )
                A.PushBack( i, i-1, -1. );
            A.PushBack( i, i, 6. );
            if( x != n-1 )
                A.PushBack( i, i+1, -1. );
            if( y != n-1 )
                A.PushBack( i, i+n, -1. );
            if( z != n-1 )
                A.PushBack( i, i+n*n, -1. );
        } 
        A.StopAssembly();
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Generating random vector x and forming y := A x...";
            std::cout.flush();
        }
        DistVector<double> x( N, comm ), y( N, comm );
        MakeUniform( x );
        MakeZeros( y );
        Multiply( 1., A, x, 0., y );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running nested dissection...";
            std::cout.flush();
        }
        const DistGraph& graph = A.Graph();
        DistSymmInfo info;
        DistSeparatorTree sepTree;
        std::vector<int> localMap;
        NestedDissection
        ( graph, localMap, sepTree, info, cutoff, numDistSeps, numSeqSeps );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

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

        if( commRank == 0 )
        {
            std::cout << "Building DistSymmFrontTree...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        DistSymmFrontTree<double> frontTree( A, localMap, sepTree, info );
        std::vector<int> inverseLocalMap;
        InvertMap( localMap, inverseLocalMap, N, comm );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running LDL^T factorization and redistribution...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        LDL( TRANSPOSE, info, frontTree );
        SetSolveMode( frontTree, FAST_2D_LDL );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Solving against y...";
            std::cout.flush();
        }
        DistNodalVector<double> yNodal;
        yNodal.Pull( inverseLocalMap, info, y );
        const double yOrigLocalNorm = elem::Nrm2( yNodal.values );
        LDLSolve( TRANSPOSE, info, frontTree, yNodal.values );
        yNodal.Push( inverseLocalMap, info, y );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
            std::cout << "Checking error in computed solution..." << std::endl;
        DistNodalVector<double> xNodal;
        xNodal.Pull( inverseLocalMap, info, x );
        const double xLocalNorm = elem::Nrm2( xNodal.values );
        const double yLocalNorm = elem::Nrm2( yNodal.values );
        elem::Axpy( -1., xNodal.values, yNodal.values );
        const double errorLocalNorm = elem::Nrm2( yNodal.values );
        const double xLocalNormSq = xLocalNorm*xLocalNorm;
        const double yLocalNormSq = yLocalNorm*yLocalNorm;
        const double yOrigLocalNormSq = yOrigLocalNorm*yOrigLocalNorm;

        const double errorLocalNormSq = errorLocalNorm*errorLocalNorm;
        double xNormSq, yNormSq, yOrigNormSq, errorNormSq;
        mpi::AllReduce( &xLocalNormSq, &xNormSq, 1, mpi::SUM, comm );
        mpi::AllReduce( &yLocalNormSq, &yNormSq, 1, mpi::SUM, comm );
        mpi::AllReduce( &yOrigLocalNormSq, &yOrigNormSq, 1, mpi::SUM, comm );
        mpi::AllReduce( &errorLocalNormSq, &errorNormSq, 1, mpi::SUM, comm );
        const double xNorm = elem::Sqrt( xNormSq );
        const double yNorm = elem::Sqrt( yNormSq );
        const double yOrigNorm = elem::Sqrt( yOrigNormSq );
        const double errorNorm = elem::Sqrt( errorNormSq );

        if( commRank == 0 )
        {
            std::cout << "|| x     ||_2 = " << xNorm << "\n"
                      << "|| xComp ||_2 = " << yNorm << "\n"
                      << "|| A x   ||_2 = " << yOrigNorm << "\n"
                      << "|| error ||_2 = " << errorNorm << std::endl;
        }
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        elem::DumpCallStack();
        cliq::DumpCallStack();
#endif
        std::ostringstream msg;
        msg << "Process " << commRank << " caught message:\n"
            << e.what() << "\n";
        std::cerr << msg.str() << std::endl;
    }

    cliq::Finalize();
    return 0;
}
