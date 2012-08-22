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
      << "SimpleMultiVectorSolve <n1> <n2> <n3> <numRhs> "
      << "[sequential=true] [numDistSeps=1] [numSeqSeps=1] [cutoff=128]\n"
      << "  n1: first dimension of n1 x n2 x n3 mesh\n"
      << "  n2: second dimension of n1 x n2 x n3 mesh\n"
      << "  n3: third dimension of n1 x n2 x n3 mesh\n"
      << "  sequential: if nonzero, then run a sequential symbolic reordering\n"
      << "  numDistSeps: number of distributed separators to try\n"
      << "  numSeqSeps: number of sequential separators to try\n"
      << "  cutoff: maximum size of leaf node\n"
      << std::endl;
}

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 5 )
    {
        if( commRank == 0 )
            Usage();
        cliq::Finalize();
        return 0;
    }
    int argNum = 1;
    const int n1 = atoi(argv[argNum++]);
    const int n2 = atoi(argv[argNum++]);
    const int n3 = atoi(argv[argNum++]);
    const int numRhs = atoi(argv[argNum++]);
    const bool sequential = ( argc>argNum ? atoi(argv[argNum++]) : true );
    const int numDistSeps = ( argc>argNum ? atoi(argv[argNum++]) : 1 );
    const int numSeqSeps = ( argc>argNum ? atoi(argv[argNum++]) : 1 );
    const int cutoff = ( argc>argNum ? atoi(argv[argNum++]) : 128 );

    try
    {
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

        if( commRank == 0 )
        {
            std::cout << "Generating random vector X and forming Y := A X...";
            std::cout.flush();
        }
        const double multiplyStart = mpi::Time();
        DistMultiVector<double> X( N, numRhs, comm ), Y( N, numRhs, comm );
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
            std::cout << "Solving...";
            std::cout.flush();
        }
        const double solveStart = mpi::Time();
        SymmetricSolve( A, Y, sequential, numDistSeps, numSeqSeps, cutoff );
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
                std::cout << "Right-hand side " << j << "\n"
                          << "------------------------------------------\n"
                          << "|| x     ||_2 = " << XNorms[j] << "\n"
                          << "|| xComp ||_2 = " << YNorms[j] << "\n"
                          << "|| A x   ||_2 = " << YOrigNorms[j] << "\n"
                          << "|| error ||_2 = " << errorNorms[j] << "\n"
                          << std::endl;
            }
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
