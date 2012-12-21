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
    typedef double R;
    typedef Complex<R> C;

    try
    {
        const int n1 = Input("--n1","first grid dimension",200);
        const int n2 = Input("--n2","second grid dimension",200);
        const double omega = Input("--omega","angular frequency",120.);
        const double damping = Input("--damping","damping parameter",7.);
        const bool analytic = Input("--analytic","analytic partitions?",true);
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
        ProcessInput();

        const int N = n1*n2;
        DistSparseMatrix<C> A( N, comm );
        C dampedOmega( omega, damping );
        const double hxInv = n1+1;
        const double hyInv = n2+1;
        const double hxInvSquared = hxInv*hxInv;
        const double hyInvSquared = hyInv*hyInv;
        const C mainTerm = 
            2*(hxInvSquared+hyInvSquared) - dampedOmega*dampedOmega;

        // Fill our portion of the 2D Helmholtz operator over the unit-square 
        // using a n1 x n2 5-point stencil in natural ordering: 
        // (x,y) at x + y*n1
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of matrix...";
            std::cout.flush();
        }
        const double fillStart = mpi::Time();
        const int firstLocalRow = A.FirstLocalRow();
        const int localHeight = A.LocalHeight();
        A.StartAssembly();
        A.Reserve( 5*localHeight );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = firstLocalRow + iLocal;
            const int x = i % n1;
            const int y = i/n1;

            A.Update( i, i, mainTerm );
            if( x != 0 )
                A.Update( i, i-1, -hxInvSquared );
            if( x != n1-1 )
                A.Update( i, i+1, -hxInvSquared );
            if( y != 0 )
                A.Update( i, i-n1, -hyInvSquared );
            if( y != n2-1 )
                A.Update( i, i+n1, -hyInvSquared );
        } 
        A.StopAssembly();
        mpi::Barrier( comm );
        const double fillStop =  mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << fillStop-fillStart << " seconds" 
                      << std::endl;

        if( print )
            A.Print("A");

        if( commRank == 0 )
        {
            std::cout << "Generating random vector x and forming y := A x...";
            std::cout.flush();
        }
        const double multiplyStart = mpi::Time();
        DistVector<C> x( N, comm ), y( N, comm );
        MakeUniform( x );
        MakeZeros( y );
        Multiply( C(1), A, x, C(0), y );
        const double yOrigNorm = Norm( y );
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
        const DistGraph& graph = A.Graph();
        DistSymmInfo info;
        DistSeparatorTree sepTree;
        DistMap map, inverseMap;
        if( analytic )
            NaturalNestedDissection
            ( n1, n2, 1, graph, map, sepTree, info, cutoff );
        else
            NestedDissection
            ( graph, map, sepTree, info, 
              sequential, numDistSeps, numSeqSeps, cutoff );
        map.FormInverse( inverseMap );
        mpi::Barrier( comm );
        const double nestedStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << nestedStop-nestedStart << " seconds"
                      << std::endl;

        const int rootSepSize = info.distNodes.back().size;
        if( commRank == 0 )
        {
            const int numDistNodes = info.distNodes.size();
            const int numLocalNodes = info.localNodes.size();
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
        const double buildStart = mpi::Time();
        DistSymmFrontTree<C> frontTree( TRANSPOSE, A, map, sepTree, info );
        mpi::Barrier( comm );
        const double buildStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << buildStop-buildStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running block LDL^T...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double ldlStart = mpi::Time();
        LDL( info, frontTree, BLOCK_LDL_2D );
        mpi::Barrier( comm );
        const double ldlStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << ldlStop-ldlStart << " seconds" 
                      << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Computing SVD of connectivity of second separator to "
                         "the root separator...";
            std::cout.flush();
        }
        const int numDistFronts = frontTree.distFronts.size();
        if( numDistFronts >= 2 && info.distNodes[numDistFronts-2].onLeft )
        {
            const double svdStart = mpi::Time();
            const DistMatrix<C>& frontL =
                frontTree.distFronts[numDistFronts-2].front2dL;
            const Grid& grid = frontL.Grid();
            const int gridRank = grid.Rank();
            const int height = frontL.Height();
            const int width = frontL.Width();
            const int minDim = std::min(height,width);
            DistMatrix<C> B( grid );
            B.LockedView( frontL, width, 0, height-width, width );
            DistMatrix<C> BCopy( B );
            DistMatrix<R,VR,STAR> singVals_VR_STAR( grid );
            elem::SingularValues( BCopy, singVals_VR_STAR );
            const R twoNorm =
                elem::Norm( singVals_VR_STAR, elem::MAX_NORM );
            DistMatrix<R,STAR,STAR> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const double svdStop = mpi::Time();
            if( gridRank == 0 )
                std::cout << "done, " << svdStop-svdStart << " seconds\n"
                          << "  two norm=" << twoNorm << "\n";
            for( double tol=1e-1; tol>=1e-10; tol/=10 )
            {
                int numRank = minDim;
                for( int j=0; j<minDim; ++j )
                {
                    if( singVals.GetLocal(j,0) <= twoNorm*tol )
                    {
                        numRank = j;
                        break;
                    }
                }
                if( gridRank == 0 )
                    std::cout << "  rank (" << tol << ")=" << numRank
                              << "/" << minDim << std::endl;
            }
        }

        if( commRank == 0 )
        {
            std::cout << "Computing SVD of the largest off-diagonal block of "
                         "numerical Green's function on root separator...";
            std::cout.flush();
        }
        {
            const double svdStart = mpi::Time();
            const DistMatrix<C>& front = frontTree.distFronts.back().front2dL;
            const Grid& grid = front.Grid();
            const int lowerHalf = rootSepSize/2;
            const int upperHalf = rootSepSize - lowerHalf;
            if( commRank == 0 )
                std::cout << "lowerHalf=" << lowerHalf 
                          << ", upperHalf=" << upperHalf << std::endl;
            DistMatrix<C> offDiagBlock;
            offDiagBlock.LockedView
            ( front, lowerHalf, 0, upperHalf, lowerHalf );
            DistMatrix<C> offDiagBlockCopy( offDiagBlock );
            DistMatrix<R,VR,STAR> singVals_VR_STAR( grid );
            elem::SingularValues( offDiagBlockCopy, singVals_VR_STAR );
            const R twoNorm = elem::Norm( singVals_VR_STAR, elem::MAX_NORM );
            const R tolerance = 1e-4;
            DistMatrix<R,STAR,STAR> singVals( singVals_VR_STAR );
            mpi::Barrier( comm );
            const double svdStop = mpi::Time();
            if( commRank == 0 ) 
                std::cout << "done, " << svdStop-svdStart << " seconds\n";
            for( double tol=1e-1; tol>=1e-10; tol/=10 )
            {
                int numRank = lowerHalf;
                for( int j=0; j<lowerHalf; ++j )
                {
                    if( singVals.GetLocal(j,0) <= twoNorm*tol )
                    {
                        numRank = j;
                        break;
                    }
                }
                if( commRank == 0 )
                    std::cout << "  rank (" << tol << ")=" << numRank
                              << "/" << lowerHalf << std::endl;
            }
        }

        if( commRank == 0 )
        {
            std::cout << "Solving against y...";
            std::cout.flush();
        }
        const double solveStart = mpi::Time();
        DistNodalVector<C> yNodal;
        yNodal.Pull( inverseMap, info, y );
        Solve( info, frontTree, yNodal.localVec );
        yNodal.Push( inverseMap, info, y );
        mpi::Barrier( comm );
        const double solveStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << solveStop-solveStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
            std::cout << "Checking error in computed solution..." << std::endl;
        const double xNorm = Norm( x );
        const double yNorm = Norm( y );
        Axpy( C(-1), x, y );
        const double errorNorm = Norm( y );
        if( commRank == 0 )
        {
            std::cout << "|| x     ||_2 = " << xNorm << "\n"
                      << "|| xComp ||_2 = " << yNorm << "\n"
                      << "|| A x   ||_2 = " << yOrigNorm << "\n"
                      << "|| error ||_2 / || x ||_2 = " 
                      << errorNorm/xNorm << "\n"
                      << "|| error ||_2 / || A x ||_2 = " 
                      << errorNorm/yOrigNorm
                      << std::endl;
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
