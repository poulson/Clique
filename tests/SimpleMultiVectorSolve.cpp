/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/
#include "clique.hpp"
using namespace cliq;

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n1 = Input("--n1","first grid dimension",30);
        const int n2 = Input("--n2","second grid dimension",30);
        const int n3 = Input("--n3","third grid dimension",30);
        const int numRhs = Input("--numRhs","number of right-hand sides",5);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const int cutoff = Input("--cutoff","cutoff for nested dissection",128);
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
    catch( ArgException& e ) { }
    catch( std::exception& e )
    {
        std::ostringstream msg;
        msg << "Process " << commRank << " caught message:\n"
            << e.what() << std::endl;
        std::cerr << msg.str();
#ifndef RELEASE
        elem::DumpCallStack();
        cliq::DumpCallStack();
#endif
    }

    cliq::Finalize();
    return 0;
}
