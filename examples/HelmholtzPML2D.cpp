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
      << "HelmholtzPML2D <nx> <ny> <omega> <Lx> <Ly> <b> <sigma> <p> "
      << "[analytic=true] [sequential=true] [cutoff=128] \n"
      << "[numDistSeps=1] [numSeqSeps=1]\n"
      << "  nx,ny: dimensions of nx x ny mesh\n"
      << "  omega: frequency of problem in radians per second\n"
      << "  Lx,Ly: dimensions of [0,Lx] x [0,Ly] domain\n"
      << "  b: width of PML in grid points\n"
      << "  sigma: coefficient for PML profile\n"
      << "  p: polynomial order for PML profile\n"
      << "  analytic: if nonzero, use an analytical reordering\n"
      << "  sequential: if nonzero, then run a sequential symbolic reordering\n"
      << "  cutoff: maximum size of leaf node\n"
      << "  numDistSeps: number of distributed separators to try\n"
      << "  numSeqSeps: number of sequential separators to try\n"
      << std::endl;
}

Complex<double> PML( double x, double w, double p, double sigma, double k )
{
#ifndef RELEASE
    if( x < 0 || x > w+1e-10 )
        throw std::logic_error("Evaluation point not in PML interval");
#endif
    const double realPart = 1.0;
    const double arg = x/w;
    const double imagPart = (sigma/w)*std::pow(arg,p)/k;
    return Complex<double>(realPart,imagPart); 
}

Complex<double> 
sInv( int j, int n, int b, double h, double p, double sigma, double k )
{
    if( j < b-1 )
        return PML( (b-1-j)*h, b*h, p, sigma, k );
    else if( j > n-b )
        return PML( (j-(n-b))*h, b*h, p, sigma, k );
    else
        return Complex<double>(1.0,0.0);
}

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    typedef double R;
    typedef Complex<R> C;

    if( argc < 9 )
    {
        if( commRank == 0 )
            Usage();
        cliq::Finalize();
        return 0;
    }
    int argNum = 1;
    const int nx = atoi(argv[argNum++]);
    const int ny = atoi(argv[argNum++]);
    const double omega = atof(argv[argNum++]);
    const double Lx = atof(argv[argNum++]);
    const double Ly = atof(argv[argNum++]);
    const int b = atoi(argv[argNum++]);
    const double sigma = atof(argv[argNum++]);
    const double p = atof(argv[argNum++]);
    const bool analytic = ( argc>argNum ? atoi(argv[argNum++]) : true );
    const bool sequential = ( argc>argNum ? atoi(argv[argNum++]) : true );
    const int cutoff = ( argc>argNum ? atoi(argv[argNum++]) : 128 );
    const int numDistSeps = ( argc>argNum ? atoi(argv[argNum++]) : 1 );
    const int numSeqSeps = ( argc>argNum ? atoi(argv[argNum++]) : 1 );

    try
    {
        const double k = omega/(2*M_PI);
        const int N = nx*ny;
        DistSparseMatrix<C> A( N, comm );
        const double hx = Lx/(nx+1);
        const double hy = Ly/(ny+1);
        const double hxSquared = hx*hx;
        const double hySquared = hy*hy;

        // Fill our portion of the 2D Helmholtz operator 
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
            const int x = i % nx;
            const int y = i/nx;

            const C s1InvL = sInv( x-1, nx, b, hx, p, sigma, k );
            const C s1InvM = sInv( x,   nx, b, hx, p, sigma, k );
            const C s1InvR = sInv( x+1, nx, b, hx, p, sigma, k );
            const C s2InvL = sInv( y-1, ny, b, hy, p, sigma, k );
            const C s2InvM = sInv( y,   ny, b, hy, p, sigma, k );
            const C s2InvR = sInv( y+1, ny, b, hy, p, sigma, k );

            const C xTop = s2InvM;
            const C xTempL = xTop/s1InvL;
            const C xTempM = xTop/s1InvM;
            const C xTempR = xTop/s1InvR;
            const C xTermL = (xTempL+xTempM) / (2*hxSquared);
            const C xTermR = (xTempM+xTempR) / (2*hxSquared);

            const C yTop = s1InvM;
            const C yTempL = yTop/s2InvL;
            const C yTempM = yTop/s2InvM;
            const C yTempR = yTop/s2InvR;
            const C yTermL = (yTempL+yTempM) / (2*hySquared);
            const C yTermR = (yTempM+yTempR) / (2*hySquared);

            const C mainTerm = (xTermL+xTermR+yTermL+yTermR)
                - omega*omega*s1InvM*s2InvM;

            A.Update( i, i, mainTerm );
            if( x != 0 )
                A.Update( i, i-1, -xTermL );
            if( x != nx-1 )
                A.Update( i, i+1, -xTermR );
            if( y != 0 )
                A.Update( i, i-nx, -yTermL );
            if( y != ny-1 )
                A.Update( i, i+nx, -yTermR );
        } 
        A.StopAssembly();
        mpi::Barrier( comm );
        const double fillStop =  mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << fillStop-fillStart << " seconds" 
                      << std::endl;

        if( commRank == 0 )
            std::cout << "Generating point-source for y..." << std::endl;
        DistVector<C> y( N, comm ), z( N, comm );
        MakeZeros( z );
        const int xSource = nx/2;
        const int ySource = ny/2;
        const int iSource = xSource + ySource*nx;
        if( iSource >= firstLocalRow && iSource < firstLocalRow+localHeight )
            z.SetLocal( iSource-firstLocalRow, Complex<double>(1.0,0.0) );
        y = z;

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
            ( nx, ny, 1, graph, map, sepTree, info, cutoff );
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

        if( mpi::CommSize( comm ) == 1 )
            y.LocalVector().Print("solution");

        if( commRank == 0 )
            std::cout << "Checking residual norm of solution..." << std::endl;
        const double bNorm = Norm( z );
        Multiply( (C)-1, A, y, (C)1, z );
        const double errorNorm = Norm( z );
        if( commRank == 0 )
        {
            std::cout << "|| b     ||_2 = " << bNorm << "\n"
                      << "|| error ||_2 / || b ||_2 = " 
                      << errorNorm/bNorm << "\n"
                      << std::endl;
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
