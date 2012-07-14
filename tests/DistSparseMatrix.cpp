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
    std::cout << "DistSparseMatrix <n>\n" 
              << "  n: size of n x n x n mesh\n"
              << "  cutoff: maximum leaf size\n"
              << std::endl;
}

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    if( argc < 3 )
    {
        if( commRank == 0 )
            Usage();
        cliq::Finalize();
        return 0;
    }
    const int n = atoi( argv[1] );
    const int cutoff = atoi( argv[2] );

    try
    {
        const int height = n*n*n;
        const int width = height;
        DistSparseMatrix<double> A( height, width, comm );

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

            double sum = 0.;
            if( z != 0 )
            {
                A.PushBack( i, i-n*n, -1. );
                sum += 1.;
            }
            if( y != 0 )
            {
                A.PushBack( i, i-n, -1. );
                sum += 1.;
            }
            if( x != 0 )
            {
                A.PushBack( i, i-1, -1. );
                sum += 1.;
            }
            if( x != n-1 )
            {
                A.PushBack( i, i+1, -1. );
                sum += 1.;
            }
            if( y != n-1 )
            {
                A.PushBack( i, i+n, -1. );
                sum += 1.;
            }
            if( z != n-1 )
            {
                A.PushBack( i, i+n*n, -1. );
                sum += 1.;
            }
            A.PushBack( i, i, sum );
        } 
        A.StopAssembly();
        mpi::Barrier( comm );
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
        NestedDissection( graph, info, sepTree, localMap, cutoff );
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
        sleep( 1 );
        DistSymmFrontTree<double> frontTree( A, localMap, sepTree, info );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;
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
