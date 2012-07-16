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
    std::cout << "PanelLDL <nx> <ny> <nz> [cutoff=16] [fact blocksize=96] "
                 "[solve blocksize=64] [useFast1d=0] [write info=0] "
                 "[basename=out]\n"
              << "<nx>: size of panel in x direction\n"
              << "<ny>: size of panel in y direction\n"
              << "<nz>: size of panel in z direction\n"
              << "[block LDL=0]: use block LDL iff != 0\n"
              << "[cutoff=16]: minimum required leaf size\n" 
              << "[fact blocksize=96]: factorization algorithmic blocksize\n"
              << "[solve blocksize=64]: solve algorithmic blocksize\n"
              << "[useFast1d=0]: use 1d distributions iff != 0\n"
              << "[write info=0]: write basic local info to file? [0/1]\n"
              << "[basename=out]: basename for each process's info file\n"
              << std::endl;
}

int ReorderedIndexRecursion
( int x, int y, int z, int nx, int ny, int nz,
  int stepsLeft, int cutoff, int offset );

int ReorderedIndex
( int x, int y, int z, int nx, int ny, int nz,
  int minBalancedDepth, int cutoff );

void CountLocalTreeSize
( int nxSub, int nySub, int nz, int cutoff, int& numNodes );

void FillElimTree
( int nx, int ny, int nz, int cutoff, mpi::Comm comm, DistSymmElimTree& eTree );

void FillDistElimTree
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, mpi::Comm comm, int distDepth, DistSymmElimTree& eTree );

void FillLocalElimTree
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int distDepth, DistSymmElimTree& eTree );

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    typedef Complex<double> F;

    if( argc < 4 )
    {
        if( commRank == 0 )        
            Usage();
        cliq::Finalize();
        return 0;
    }

    int argNum = 1;
    const int nx = atoi(argv[argNum++]);
    const int ny = atoi(argv[argNum++]);
    const int nz = atoi(argv[argNum++]);
    const bool useBlockLDL = ( argc >= 5 ? atoi(argv[argNum++]) : false );
    const int cutoff = ( argc >= 6 ? atoi(argv[argNum++]) : 16 );
    const int factBlocksize = ( argc >= 7 ? atoi(argv[argNum++]) : 96 );
    const int solveBlocksize = ( argc >= 8 ? atoi(argv[argNum++]) : 64 );
    const bool useFast1d = ( argc >= 9 ? atoi(argv[argNum++]) : 0 );
    const bool writeInfo = ( argc >= 10 ? atoi(argv[argNum++]) : 0 );
    const char* basename = ( argc >= 11 ? argv[argNum++] : "out" );
    if( commRank == 0 )
        std::cout << "(nx,ny,nz)=(" << nx << "," << ny << "," << nz << ")\n"
                  << "cutoff=" << cutoff << "\n"
                  << "factBlocksize=" << factBlocksize << "\n"
                  << "solveBlocksize=" << solveBlocksize << std::endl;

    try
    {
        std::ofstream infoFile;
        if( writeInfo )
        {
            std::ostringstream filename;
            filename << basename << "-" << commRank << ".dat";
            infoFile.open( filename.str().c_str() );
        }

        // Fill the distributed portion of the original structure
        mpi::Barrier( comm );
        const double initStartTime = mpi::Time();
        DistSymmElimTree eTree;
        FillElimTree( nx, ny, nz, cutoff, comm, eTree );
        const int numLocalNodes = eTree.localNodes.size();
        const int numDistNodes = eTree.distNodes.size();
        mpi::Barrier( comm );
        const double initStopTime = mpi::Time();
        if( commRank == 0 )
            std::cout << "Init time: " << initStopTime-initStartTime << " secs"
                      << std::endl;

        if( writeInfo )
        {
            infoFile << "Local original structure sizes:\n";
            for( int s=0; s<numLocalNodes; ++s )
                infoFile << "  " << s << ": node=" 
                         << eTree.localNodes[s]->size << ", lower="
                         << eTree.localNodes[s]->lowerStruct.size() << "\n"; 
            infoFile << "\nDist original structure sizes:\n";
            for( int s=0; s<numDistNodes; ++s )
                infoFile << "  " << s << ": node="
                         << eTree.distNodes[s].size << ", lower="
                         << eTree.distNodes[s].lowerStruct.size() << "\n";
            infoFile << std::endl;
        }

        // Call the symbolic factorization routine
        const double analysisStartTime = mpi::Time();
        DistSymmInfo info;
        SymmetricAnalysis( eTree, info, true );
        mpi::Barrier( comm );
        const double analysisStopTime = mpi::Time();
        if( commRank == 0 )
            std::cout << "Analysis time: " 
                      << analysisStopTime-analysisStartTime << " secs"
                      << std::endl;

        if( writeInfo )
        {
            infoFile << "Local factor structure sizes:\n";
            for( int s=0; s<numLocalNodes; ++s )
                infoFile << "  " << s << ": node=" 
                         << info.localNodes[s].size << ", "
                         << " lower=" 
                         << info.localNodes[s].lowerStruct.size() << "\n";
            infoFile << "\n"
                     << "Dist factor structure sizes:\n";
            for( int s=0; s<numDistNodes; ++s )
                infoFile << "  " << s << ": node="
                         << info.distNodes[s].size << ", "
                         << " lower="
                         << info.distNodes[s].lowerStruct.size() << "\n";
            infoFile << std::endl;

            infoFile << "\n"
                     << "Local height 1d: " 
                     << info.distNodes.back().localOffset1d + 
                        info.distNodes.back().localSize1d << std::endl;
        }

        const int numPanels = 10;
        std::vector<DistSymmFrontTree<F>*> symmFrontTrees;
        for( int panel=0; panel<numPanels; ++panel )
        {
            if( commRank == 0 )
                std::cout << "Panel " << panel << " of " << numPanels 
                          << std::endl;
            if( writeInfo )
                infoFile << "Panel " << panel << " of " << numPanels 
                         << std::endl;
            // Directly initialize the frontal matrices with the original 
            // sparse matrix (for now, use an original matrix equal to identity)
            const double fillStartTime = mpi::Time();
            symmFrontTrees.push_back( new DistSymmFrontTree<F> );
            DistSymmFrontTree<F>& L = *symmFrontTrees.back();
            L.localFronts.resize( info.localNodes.size() );
            for( int s=0; s<numLocalNodes; ++s )
            {
                const LocalSymmNodeInfo& node = info.localNodes[s];
                Matrix<F>& frontL = L.localFronts[s].frontL;

                const int frontSize = node.size+node.lowerStruct.size();
                elem::Uniform( frontSize, node.size, frontL );
                if( writeInfo )
                    frontL.Print( infoFile, "frontL local" );
            }
            L.mode = NORMAL_2D;
            L.distFronts.resize( numDistNodes );
            InitializeDistLeaf( info, L );
            for( int s=1; s<numDistNodes; ++s )
            {
                const DistSymmNodeInfo& node = info.distNodes[s];
                DistMatrix<F>& front2dL = L.distFronts[s].front2dL;

                front2dL.SetGrid( *node.grid );
                const int frontSize = node.size+node.lowerStruct.size();
                front2dL.Align( 0, 0 );
                elem::Uniform( frontSize, node.size, front2dL );
                if( writeInfo )
                    front2dL.Print( infoFile, "frontL dist" );
            }
            mpi::Barrier( comm );
            const double fillStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Fill time: " << fillStopTime-fillStartTime
                          << " secs" << std::endl;

            // Generate a random RHS, multiply it by our matrix, and then 
            // compare the original RHS against the solution.
            const int NUM_RHS = 1;
            const double makeRhsStartTime = mpi::Time();
            const int localHeight1d = 
                info.distNodes.back().localOffset1d + 
                info.distNodes.back().localSize1d;
            Matrix<F> localX;
            elem::Uniform( localHeight1d, NUM_RHS, localX );
            Matrix<F> localYLower = localX;
            SetSolveMode( L, NORMAL_1D );
            LowerMultiply( NORMAL, NON_UNIT, -1, info, L, localYLower );
            Matrix<F> localY = localX;
            LowerMultiply( TRANSPOSE, NON_UNIT, 0, info, L, localY );
            elem::Axpy( (F)1, localYLower, localY );
            localYLower.Empty();
            SetSolveMode( L, NORMAL_2D );
            mpi::Barrier( comm );
            const double makeRhsStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Make RHS time: " 
                          << makeRhsStopTime-makeRhsStartTime
                          << " secs" << std::endl;
            if( writeInfo )
            {
                localX.Print( infoFile, "localX" );
                localY.Print( infoFile, "localY" );
            }
            const double myYNorm = elem::Norm( localY );
            double YNorm;
            mpi::Reduce( &myYNorm, &YNorm, 1, mpi::SUM, 0, comm );

            // Call the numerical factorization routine
            elem::SetBlocksize( factBlocksize );
            const double factStartTime = mpi::Time();
            if( useBlockLDL )
                BlockLDL( TRANSPOSE, info, L );
            else
                LDL( TRANSPOSE, info, L );
            mpi::Barrier( comm );
            const double factStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Factorization time: " 
                          << factStopTime-factStartTime
                          << " secs" << std::endl;

            if( !useBlockLDL )
            {
                // Invert the diagonal blocks for faster solves
                const double redistStartTime = mpi::Time();
                if( useFast1d )
                    SetSolveMode( L, FAST_1D_LDL );
                else
                    SetSolveMode( L, FAST_2D_LDL );
                mpi::Barrier( comm );
                const double redistStopTime = mpi::Time();
                if( commRank == 0 )
                    std::cout << "Redistribution time: " 
                              << redistStopTime-redistStartTime 
                              << " secs" << std::endl;
            }

            // Solve
            elem::SetBlocksize( solveBlocksize );
            mpi::Barrier( comm );
            const double solveStartTime = mpi::Time();
            if( useBlockLDL )
                BlockLDLSolve( TRANSPOSE, info, L, localY );
            else
                LDLSolve( TRANSPOSE, info, L, localY );
            mpi::Barrier( comm );
            const double solveStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Solve time: " << solveStopTime-solveStartTime
                          << " secs" << std::endl;

            if( writeInfo )
                localY.Print( infoFile, "localY final" );
            elem::Axpy( (F)-1, localX, localY );
            const double myErrorNorm = elem::Norm( localY );
            double errorNorm;
            mpi::Reduce( &myErrorNorm, &errorNorm, 1, mpi::SUM, 0, comm );
            if( commRank == 0 )
            {
                std::cout << "||y||_2: " << YNorm << "\n"
                          << "||e||_2: " << errorNorm << "\n" 
                          << "||e||_2/||y||_2: " << errorNorm/YNorm << "\n"
                          << std::endl;
            }
            if( writeInfo )
                localY.Print( infoFile, "error final" );
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
 
int ReorderedIndexRecursion
( int x, int y, int z, int nx, int ny, int nz,
  int stepsLeft, int cutoff, int offset )
{
    const int size = nx*ny*nz;
#ifndef RELEASE
    if( stepsLeft != 0 && size == 0 )
        throw std::logic_error("Null node in the upper tree");
#endif
    if( stepsLeft == 0 && size <= cutoff )
    {
        // We have satisfied the nested dissection constraints
        return offset + (x+y*nx+z*nx*ny);
    }
    else if( nx >= ny )
    {
        // Partition the X dimension
        const int middle = (nx-1)/2;
        if( x < middle )
        {
            return ReorderedIndexRecursion
            ( x, y, z, middle, ny, nz, 
              std::max(stepsLeft-1,0), cutoff, offset );
        }
        else if( x == middle )
        {
            return offset + std::max(nx-1,0)*ny*nz + (y+z*ny);
        }
        else // x > middle
        {
            return ReorderedIndexRecursion
            ( x-middle-1, y, z, std::max(nx-middle-1,0), ny, nz,
              std::max(stepsLeft-1,0), cutoff, offset+middle*ny*nz );
        }
    }
    else
    {
        // Partition the Y dimension
        const int middle = (ny-1)/2;
        if( y < middle )
        {
            return ReorderedIndexRecursion
            ( x, y, z, nx, middle, nz,
              std::max(stepsLeft-1,0), cutoff, offset );
        }
        else if( y == middle )
        {
            return offset + nx*std::max(ny-1,0)*nz + (x+z*nx);
        }
        else // y > middle 
        {
            return ReorderedIndexRecursion
            ( x, y-middle-1, z, nx, std::max(ny-middle-1,0), nz,
              std::max(stepsLeft-1,0), cutoff, offset+nx*middle*nz );
        }
    }
}

int ReorderedIndex
( int x, int y, int z, int nx, int ny, int nz,
  int minBalancedDepth, int cutoff )
{
#ifndef RELEASE    
    cliq::PushCallStack("ReorderedIndex");
#endif
    int index = ReorderedIndexRecursion
    ( x, y, z, nx, ny, nz, minBalancedDepth, cutoff, 0 );
#ifndef RELEASE
    cliq::PopCallStack();
#endif
    return index;
}

void FillElimTree
( int nx, int ny, int nz, int cutoff, mpi::Comm comm, 
  DistSymmElimTree& eTree )
{
#ifndef RELEASE
    cliq::PushCallStack("FillElimTree");
#endif
    const int distDepth = DistributedDepth( comm );
    int nxSub=nx, nySub=ny, xOffset=0, yOffset=0;
    FillDistElimTree
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      comm, distDepth, eTree );
    FillLocalElimTree
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      distDepth, eTree );
#ifndef RELEASE
    cliq::PopCallStack();
#endif
}
  
void FillDistElimTree
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, mpi::Comm comm, int distDepth, DistSymmElimTree& eTree )
{
#ifndef RELEASE
    cliq::PushCallStack("FillDistElimTree");
#endif
    const int numDistNodes = distDepth+1;

    eTree.distNodes.resize( numDistNodes );
    mpi::CommDup( comm, eTree.distNodes.back().comm );

    // Fill the distributed nodes
    for( int s=numDistNodes-1; s>0; --s )
    {
        DistSymmNode& node = eTree.distNodes[s];
        DistSymmNode& childNode = eTree.distNodes[s-1];

        const int nodeCommRank = mpi::CommRank( node.comm );
        const int nodeCommSize = mpi::CommSize( node.comm );
        const int leftTeamSize = nodeCommSize/2;

        const bool onLeft = ( nodeCommRank < leftTeamSize );
        const int childNodeCommRank = 
            ( onLeft ? nodeCommRank : nodeCommRank-leftTeamSize );
        mpi::CommSplit( node.comm, onLeft, childNodeCommRank, childNode.comm );
        childNode.onLeft = onLeft;

        if( nxSub >= nySub )
        {
            // Form the structure of a partition of the X dimension
            const int middle = (nxSub-1)/2;
            node.size = nySub*nz;
            node.offset = 
                ReorderedIndex
                ( xOffset+middle, yOffset, 0, nx, ny, nz, 
                  distDepth, cutoff );

            // Allocate space for the lower structure
            int numJoins = 0;
            if( yOffset-1 >= 0 )
                ++numJoins;
            if( yOffset+nySub < ny )
                ++numJoins;
            node.lowerStruct.resize( numJoins*nz );

            // Fill the (unsorted) lower structure
            int joinOffset = 0;
            if( yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    node.lowerStruct[i] = ReorderedIndex
                    ( xOffset+middle, yOffset-1, i, nx, ny, nz, 
                      distDepth, cutoff );
                joinOffset += nz;
            }
            if( yOffset+nySub < ny )
            {
                for( int i=0; i<nz; ++i )
                    node.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+middle, yOffset+nySub, i, nx, ny, nz, 
                      distDepth, cutoff );
            }

            // Sort the lower structure
            std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );
            
            // Pick the new offsets and sizes based upon our rank
            if( onLeft )
            {
                xOffset = xOffset;
                nxSub = middle;
            }
            else
            {
                xOffset = xOffset+middle+1;
                nxSub = std::max(nxSub-middle-1,0);
            }
        }
        else
        {
            // Form the structure of a partition of the Y dimension
            const int middle = (nySub-1)/2;
            node.size = nxSub*nz;
            node.offset = 
                ReorderedIndex
                ( xOffset, yOffset+middle, 0, nx, ny, nz, 
                  distDepth, cutoff );

            // Allocate space for the lower structure
            int numJoins = 0;
            if( xOffset-1 >= 0 )
                ++numJoins;
            if( xOffset+nxSub < nx )
                ++numJoins;
            node.lowerStruct.resize( numJoins*nz );

            // Fill the (unsorted) lower structure
            int joinOffset = 0;
            if( xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    node.lowerStruct[i] = ReorderedIndex
                    ( xOffset-1, yOffset+middle, i, nx, ny, nz, 
                      distDepth, cutoff );
                joinOffset += nz;
            }
            if( xOffset+nxSub < nx )
            {
                for( int i=0; i<nz; ++i )
                    node.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+nxSub, yOffset+middle, i, nx, ny, nz,
                      distDepth, cutoff );
            }

            // Sort the lower structure
            std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );

            // Pick the new offsets and sizes based upon our rank
            if( onLeft )
            {
                yOffset = yOffset;
                nySub = middle;
            }
            else
            {
                yOffset = yOffset+middle+1;
                nySub = std::max(nySub-middle-1,0);
            }
        }
    }

    // Fill the bottom node, which is only owned by a single process
    DistSymmNode& node = eTree.distNodes[0];
    if( nxSub*nySub*nz <= cutoff )
    {
        node.size = nxSub*nySub*nz;
        node.offset = ReorderedIndex
            ( xOffset, yOffset, 0, nx, ny, nz, distDepth, cutoff );

        // Count, allocate, and fill the lower struct
        int joinSize = 0;
        if( xOffset-1 >= 0 )
            joinSize += nySub*nz;
        if( xOffset+nxSub < nx )
            joinSize += nySub*nz;
        if( yOffset-1 >= 0 )
            joinSize += nxSub*nz;
        if( yOffset+nySub < ny )
            joinSize += nxSub*nz;
        node.lowerStruct.resize( joinSize );

        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nySub; ++j )
                    node.lowerStruct[i*nySub+j] = ReorderedIndex
                    ( xOffset-1, yOffset+j, i, 
                      nx, ny, nz, distDepth, cutoff );
            joinOffset += nySub*nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nySub; ++j )
                    node.lowerStruct[joinOffset+i*nySub+j] = ReorderedIndex
                    ( xOffset+nxSub, yOffset+j, i, 
                      nx, ny, nz, distDepth, cutoff );
            joinOffset += nySub*nz;
        }
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nxSub; ++j )
                    node.lowerStruct[joinOffset+i*nxSub+j] = ReorderedIndex
                    ( xOffset+j, yOffset-1, i, 
                      nx, ny, nz, distDepth, cutoff );
            joinOffset += nxSub*nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                for( int j=0; j<nxSub; ++j )
                    node.lowerStruct[joinOffset+i*nxSub+j] = ReorderedIndex
                    ( xOffset+j, yOffset+nySub, i, 
                      nx, ny, nz, distDepth, cutoff );
        }

        // Sort the lower structure
        std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );
    }
    else if( nxSub >= nySub )
    {
        // Form the structure of a partition of the X dimension
        const int middle = (nxSub-1)/2;
        node.size = nySub*nz;
        node.offset = 
            ReorderedIndex
            ( xOffset+middle, yOffset, 0, nx, ny, nz, distDepth, cutoff );

        // Allocate space for the lower structure
        int numJoins = 0;
        if( yOffset-1 >= 0 )
            ++numJoins;
        if( yOffset+nySub < ny )
            ++numJoins;
        node.lowerStruct.resize( numJoins*nz );

        // Fill the (unsorted) lower structure
        int joinOffset = 0;
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                node.lowerStruct[i] = ReorderedIndex
                ( xOffset+middle, yOffset-1, i, nx, ny, nz, 
                  distDepth, cutoff );
            joinOffset += nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                node.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+middle, yOffset+nySub, i, nx, ny, nz, 
                  distDepth, cutoff );
        }

        // Sort the lower structure
        std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );
    }
    else
    {
        // Form the structure of a partition of the Y dimension
        const int middle = (nySub-1)/2;
        node.size = nxSub*nz;
        node.offset = 
            ReorderedIndex
            ( xOffset, yOffset+middle, 0, nx, ny, nz, distDepth, cutoff );

        // Allocate space for the lower structure
        int numJoins = 0;
        if( xOffset-1 >= 0 )
            ++numJoins;
        if( xOffset+nxSub < nx )
            ++numJoins;
        node.lowerStruct.resize( numJoins*nz );

        // Fill the (unsorted) lower structure
        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                node.lowerStruct[i] = ReorderedIndex
                ( xOffset-1, yOffset+middle, i, nx, ny, nz, 
                  distDepth, cutoff );
            joinOffset += nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                node.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+nxSub, yOffset+middle, i, nx, ny, nz,
                  distDepth, cutoff );
        }

        // Sort the lower structure
        std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );
    }
#ifndef RELEASE
    cliq::PopCallStack();
#endif
}

// Even though this is only used within the following function, it must be
// declared outside of the function in order to be used as a template argument
// (for std::stack) since template arguments must have external linkage.
//
// This data holds the bounding box information.
struct Box
{
    int parentIndex, nx, ny, xOffset, yOffset;
    bool leftChild;
};

void FillLocalElimTree
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int distDepth, DistSymmElimTree& eTree )
{
#ifndef RELEASE
    cliq::PushCallStack("FillLocalElimTree");
#endif
    // First count the depth, resize, and then run the actual fill
    int numNodes = 0;
    CountLocalTreeSize( nxSub, nySub, nz, cutoff, numNodes );
    eTree.localNodes.resize( numNodes );
    
    // Initialize with the root's box
    std::stack<Box> boxStack;
    {
        Box box;
        box.parentIndex = -1;
        box.nx = nxSub;
        box.ny = nySub;
        box.xOffset = xOffset;
        box.yOffset = yOffset;
        box.leftChild = false;
        boxStack.push(box);
    }

    // Fill the local tree
    for( int s=numNodes-1; s>=0; --s )
    {
        Box box = boxStack.top();
        boxStack.pop();

        eTree.localNodes[s] = new LocalSymmNode;
        LocalSymmNode& node = *eTree.localNodes[s];
        node.parent = box.parentIndex;
        if( node.parent != -1 )
        {
            if( box.leftChild )
                eTree.localNodes[node.parent]->children[0] = s;
            else
                eTree.localNodes[node.parent]->children[1] = s;
        }

        if( box.nx*box.ny*nz <= cutoff )
        {
            node.size = box.nx*box.ny*nz;
            node.offset = ReorderedIndex
                ( box.xOffset, box.yOffset, 0, nx, ny, nz, 
                  distDepth, cutoff );
            node.children.clear();

            // Count, allocate, and fill the lower struct
            int joinSize = 0;
            if( box.xOffset-1 >= 0 )
                joinSize += box.ny*nz;
            if( box.xOffset+box.nx < nx )
                joinSize += box.ny*nz;
            if( box.yOffset-1 >= 0 )
                joinSize += box.nx*nz;
            if( box.yOffset+box.ny < ny )
                joinSize += box.nx*nz;
            node.lowerStruct.resize( joinSize );

            int joinOffset = 0;
            if( box.xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.ny; ++j )
                        node.lowerStruct[i*box.ny+j] = ReorderedIndex
                        ( box.xOffset-1, box.yOffset+j, i, 
                          nx, ny, nz, distDepth, cutoff );
                joinOffset += box.ny*nz;
            }
            if( box.xOffset+box.nx < nx )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.ny; ++j )
                        node.lowerStruct[joinOffset+i*box.ny+j] = ReorderedIndex
                        ( box.xOffset+box.nx, box.yOffset+j, i, 
                          nx, ny, nz, distDepth, cutoff );
                joinOffset += box.ny*nz;
            }
            if( box.yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.nx; ++j )
                        node.lowerStruct[joinOffset+i*box.nx+j] = ReorderedIndex
                        ( box.xOffset+j, box.yOffset-1, i, 
                          nx, ny, nz, distDepth, cutoff );
                joinOffset += box.nx*nz;
            }
            if( box.yOffset+box.ny < ny )
            {
                for( int i=0; i<nz; ++i )
                    for( int j=0; j<box.nx; ++j )
                        node.lowerStruct[joinOffset+i*box.nx+j] = ReorderedIndex
                        ( box.xOffset+j, box.yOffset+box.ny, i,
                          nx, ny, nz, distDepth, cutoff );
            }

            // Sort the lower structure
            std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );
        }
        else
        {
            node.children.resize(2);
            if( box.nx >= box.ny )
            {
                // Partition the X dimension (this is the separator)
                const int middle = (box.nx-1)/2;
                node.size = box.ny*nz;
                node.offset = ReorderedIndex
                    ( box.xOffset+middle, box.yOffset, 0, nx, ny, nz,
                      distDepth, cutoff );

                // Count, allocate, and fill the lower struct
                int numJoins = 0;
                if( box.yOffset-1 >= 0 )
                    ++numJoins;
                if( box.yOffset+box.ny < ny )
                    ++numJoins;
                node.lowerStruct.resize( numJoins*nz );

                int joinOffset = 0;
                if( box.yOffset-1 >= 0 )
                {
                    for( int i=0; i<nz; ++i )
                        node.lowerStruct[i] = ReorderedIndex
                        ( box.xOffset+middle, box.yOffset-1, i, nx, ny, nz,
                          distDepth, cutoff );
                    joinOffset += nz;
                }
                if( box.yOffset+box.ny < ny )
                {
                    for( int i=0; i<nz; ++i )
                        node.lowerStruct[joinOffset+i] = ReorderedIndex
                        ( box.xOffset+middle, box.yOffset+box.ny, i, nx, ny, nz,
                          distDepth, cutoff );
                }

                // Sort the lower structure
                std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );

                // Push the left child box onto the stack
                Box leftBox;
                leftBox.parentIndex = s;
                leftBox.nx = middle;
                leftBox.ny = box.ny;
                leftBox.xOffset = box.xOffset;
                leftBox.yOffset = box.yOffset;
                leftBox.leftChild = true;
                boxStack.push( leftBox );

                // Push the right child box onto the stack
                Box rightBox;
                rightBox.parentIndex = s;
                rightBox.nx = std::max(box.nx-middle-1,0);
                rightBox.ny = box.ny;
                rightBox.xOffset = box.xOffset+middle+1;
                rightBox.yOffset = box.yOffset;
                rightBox.leftChild = false;
                boxStack.push( rightBox );
            }
            else 
            {
                // Partition the Y dimension (this is the separator)
                const int middle = (box.ny-1)/2;
                node.size = box.nx*nz;
                node.offset = ReorderedIndex
                    ( box.xOffset, box.yOffset+middle, 0, nx, ny, nz,
                      distDepth, cutoff );

                // Count, allocate, and fill the lower struct
                int numJoins = 0;
                if( box.xOffset-1 >= 0 )
                    ++numJoins;
                if( box.xOffset+box.nx < nx )
                    ++numJoins;
                node.lowerStruct.resize( numJoins*nz );

                int joinOffset = 0;
                if( box.xOffset-1 >= 0 )
                {
                    for( int i=0; i<nz; ++i )
                        node.lowerStruct[i] = ReorderedIndex
                        ( box.xOffset-1, box.yOffset+middle, i, nx, ny, nz,
                          distDepth, cutoff );
                    joinOffset += nz;
                }
                if( box.xOffset+box.nx < nx )
                {
                    for( int i=0; i<nz; ++i )
                        node.lowerStruct[joinOffset+i] = ReorderedIndex
                        ( box.xOffset+box.nx, box.yOffset+middle, i, nx, ny, nz,
                          distDepth, cutoff );
                }

                // Sort the lower structure
                std::sort( node.lowerStruct.begin(), node.lowerStruct.end() );

                // Push the left child box onto the stack
                Box leftBox;
                leftBox.parentIndex = s;
                leftBox.nx = box.nx;
                leftBox.ny = middle;
                leftBox.xOffset = box.xOffset;
                leftBox.yOffset = box.yOffset;
                leftBox.leftChild = true;
                boxStack.push( leftBox );
                
                // Push the right child box onto the stack
                Box rightBox;
                rightBox.parentIndex = s;
                rightBox.nx = box.nx;
                rightBox.ny = std::max(box.ny-middle-1,0);
                rightBox.xOffset = box.xOffset;
                rightBox.yOffset = box.yOffset+middle+1;
                rightBox.leftChild = false;
                boxStack.push( rightBox );
            }
        }
    }
#ifndef RELEASE
    cliq::PopCallStack();
#endif
}

void CountLocalTreeSize
( int nx, int ny, int nz, int cutoff, int& numNodes )
{
    ++numNodes;
    const int size = nx*ny*nz;
    if( size <= cutoff )
    {
        // no-op
    }
    else if( nx >= ny )
    {
        // Partition the X dimension
        const int middle = (nx-1)/2;

        // Recurse on the left partition
        CountLocalTreeSize( middle, ny, nz, cutoff, numNodes );

        // Recurse on the right partition
        CountLocalTreeSize( std::max(nx-middle-1,0), ny, nz, cutoff, numNodes );
    }
    else
    {
        // Partition the Y dimension
        const int middle = (ny-1)/2;

        // Recurse on the top partition
        CountLocalTreeSize( nx, middle, nz, cutoff, numNodes );

        // Recurse on the bottom partition
        CountLocalTreeSize( nx, std::max(ny-middle-1,0), nz, cutoff, numNodes );
    }
}
