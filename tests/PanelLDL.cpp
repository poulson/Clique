/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011 Jack Poulson, Lexing Ying, and 
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
using namespace elemental;

void Usage()
{
    std::cout << "PanelLDL <nx> <ny> <nz> <cutoff>\n"
              << "<nx>: size of panel in x direction\n"
              << "<ny>: size of panel in y direction\n"
              << "<nz>: size of panel in z direction\n"
              << "<cutoff>: minimum required leaf size\n" << std::endl;
}

int ReorderedIndexRecursion
( int x, int y, int z, int nx, int ny, int nz,
  int stepsLeft, int cutoff, int offset );

int ReorderedIndex
( int x, int y, int z, int nx, int ny, int nz,
  int minBalancedDepth, int cutoff );

void CountLocalTreeSize
( int nxSub, int nySub, int nz, int cutoff, int& numSupernodes );

void FillOrigStruct
( int nx, int ny, int nz, int cutoff, mpi::Comm comm, int log2CommSize,
  clique::symbolic::LocalSymmOrig& localOrig, 
  clique::symbolic::DistSymmOrig& distOrig );

void FillDistOrigStruct
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, mpi::Comm comm, int log2CommSize,
  clique::symbolic::DistSymmOrig& SOrig );

void FillLocalOrigStruct
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int log2CommSize,
  clique::symbolic::LocalSymmOrig& SOrig );

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );
    typedef std::complex<double> F;

    // Ensure that we have a power of two number of processes
    unsigned temp = commSize;
    unsigned log2CommSize = 0;
    while( temp >>= 1 )
        ++log2CommSize;
    if( 1u<<log2CommSize != commSize )
    {
        if( commRank == 0 )
        {
            std::cerr << "Must use a power of two number of processes" 
                      << std::endl;
        }
        clique::Finalize();
        return 0;
    }

    if( argc < 5 )
    {
        if( commRank == 0 )        
            Usage();
        clique::Finalize();
        return 0;
    }

    int argNum = 1;
    const int nx = atoi( argv[argNum++] );
    const int ny = atoi( argv[argNum++] );
    const int nz = atoi( argv[argNum++] );
    const int cutoff = atoi( argv[argNum++] );
    if( commRank == 0 )
        std::cout << "(nx,ny,nz)=(" << nx << "," << ny << "," << nz << ")\n"
                  << "cutoff=" << cutoff << std::endl;

    try
    {
        // Fill the distributed portion of the original structure
        clique::symbolic::LocalSymmOrig localSOrig;
        clique::symbolic::DistSymmOrig distSOrig;
        FillOrigStruct
        ( nx, ny, nz, cutoff, comm, log2CommSize, localSOrig, distSOrig );

        // Call the symbolic factorization routine
        clique::symbolic::DistSymmFact distS;
        clique::symbolic::LocalSymmFact localS;
        clique::symbolic::SymmetricFactorization
        ( localSOrig, distSOrig, localS, distS, true );

        // Directly initialize the frontal matrices with the original 
        // sparse matrix (for now, use an original matrix equal to identity)
        clique::numeric::LocalSymmFact<F> localL;
        localL.supernodes.resize( localS.supernodes.size() );
        for( int s=0; s<localS.supernodes.size(); ++s )
        {
            const clique::symbolic::LocalSymmFactSupernode& symbSN  =
                localS.supernodes[s];
            clique::numeric::LocalSymmFactSupernode<F>& sn = 
                localL.supernodes[s];

            const int frontSize = symbSN.size+symbSN.lowerStruct.size();
            sn.front.ResizeTo( frontSize, frontSize );
            sn.front.SetToZero();
            Matrix<F> frontTL;
            frontTL.View( sn.front, 0, 0, symbSN.size, symbSN.size );
            frontTL.SetToIdentity();
        }
        clique::numeric::DistSymmFact<F> distL;
        distL.mode = clique::MANY_RHS;
        distL.supernodes.resize( log2CommSize+1 );
        for( int s=0; s<log2CommSize+1; ++s )
        {
            const clique::symbolic::DistSymmFactSupernode& symbSN = 
                distS.supernodes[s];
            clique::numeric::DistSymmFactSupernode<F>& sn = distL.supernodes[s];

            sn.front2d.SetGrid( *symbSN.grid );
            const int frontSize = symbSN.size+symbSN.lowerStruct.size();
            sn.front2d.ResizeTo( frontSize, frontSize );
            sn.front2d.SetToZero();
            DistMatrix<F,MC,MR> frontTL;
            frontTL.View( sn.front2d, 0, 0, symbSN.size, symbSN.size );
            frontTL.SetToIdentity();
        }

        // Call the numerical factorization routine
        clique::numeric::LDL( ADJOINT, localS, distS, localL, distL );

        // Set up the properly ordered RHS and call a solve routine
        Matrix<F> localX;
        // TODO
        //clique::numeric::LDLSolve
        //( ADJOINT, localS, distS, localL, distL, localX, true );
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        elemental::DumpCallStack();
        clique::DumpCallStack();
#endif
        std::ostringstream msg;
        msg << "Process " << commRank << " caught message:\n"
            << e.what() << "\n";
        std::cerr << msg.str() << std::endl;
    }

    clique::Finalize();
    return 0;
}
 
int ReorderedIndexRecursion
( int x, int y, int z, int nx, int ny, int nz,
  int stepsLeft, int cutoff, int offset )
{
    const int size = nx*ny*nz;
#ifndef RELEASE
    if( stepsLeft != 0 && size == 0 )
        throw std::logic_error("Null supernode in the upper tree");
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
    clique::PushCallStack("ReorderedIndex");
#endif
    int index = ReorderedIndexRecursion
    ( x, y, z, nx, ny, nz, minBalancedDepth, cutoff, 0 );
#ifndef RELEASE
    clique::PopCallStack();
#endif
    return index;
}

void FillOrigStruct
( int nx, int ny, int nz, int cutoff, mpi::Comm comm, int log2CommSize,
  clique::symbolic::LocalSymmOrig& localSOrig, 
  clique::symbolic::DistSymmOrig& distSOrig )
{
#ifndef RELEASE
    clique::PushCallStack("FillOrigStruct");
#endif
    int nxSub=nx, nySub=ny, xOffset=0, yOffset=0;
    FillDistOrigStruct
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      comm, log2CommSize, distSOrig );
    FillLocalOrigStruct
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      log2CommSize, localSOrig );
#ifndef RELEASE
    clique::PopCallStack();
#endif
}
  
void FillDistOrigStruct
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, mpi::Comm comm, int log2CommSize, 
  clique::symbolic::DistSymmOrig& SOrig )
{
#ifndef RELEASE
    clique::PushCallStack("FillDistOrigStruct");
#endif
    const int commRank = mpi::CommRank( comm );
    SOrig.comm = comm;
    SOrig.supernodes.resize( log2CommSize+1 );
    // Fill the distributed nodes
    for( int s=log2CommSize; s>0; --s )
    {
        clique::symbolic::DistSymmOrigSupernode& sn = SOrig.supernodes[s];
        const int powerOfTwo = 1u<<(s-1);
        const bool onLeft = (commRank&powerOfTwo) == 0;
        if( nxSub >= nySub )
        {
            // Form the structure of a partition of the X dimension
            const int middle = (nxSub-1)/2;
            sn.size = nySub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset+middle, yOffset, 0, nx, ny, nz, 
                  log2CommSize, cutoff );

            // Allocate space for the lower structure
            int numJoins = 0;
            if( yOffset-1 >= 0 )
                ++numJoins;
            if( yOffset+nySub < ny )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            // Fill the (unsorted) lower structure
            int joinOffset = 0;
            if( yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( xOffset+middle, yOffset-1, i, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( yOffset+nySub < ny )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+middle, yOffset+nySub, i, nx, ny, nz, 
                      log2CommSize, cutoff );
            }

            // Sort the lower structure
            std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
            
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
            sn.size = nxSub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset, yOffset+middle, 0, nx, ny, nz, 
                  log2CommSize, cutoff );

            // Allocate space for the lower structure
            int numJoins = 0;
            if( xOffset-1 >= 0 )
                ++numJoins;
            if( xOffset+nxSub < nx )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            // Fill the (unsorted) lower structure
            int joinOffset = 0;
            if( xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( xOffset-1, yOffset+middle, i, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( xOffset+nxSub < nx )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+nxSub, yOffset+middle, i, nx, ny, nz,
                      log2CommSize, cutoff );
            }

            // Sort the lower structure
            std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );

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
    clique::symbolic::DistSymmOrigSupernode& sn = SOrig.supernodes[0];
    if( nxSub*nySub*nz <= cutoff )
    {
        sn.size = nxSub*nySub*nz;
        sn.offset = ReorderedIndex
            ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );

        // Count, allocate, and fill the lower struct
        int numJoins = 0;
        if( xOffset-1 >= 0 )
            ++numJoins;
        if( xOffset+nxSub < nx )
            ++numJoins;
        if( yOffset-1 >= 0 )
            ++numJoins;
        if( yOffset+nySub < ny )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset-1, yOffset, i, nx, ny, nz, log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+nxSub, yOffset, i, nx, ny, nz, log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset, yOffset-1, i, nx, ny, nz, log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset, yOffset+nySub, i, nx, ny, nz, log2CommSize, cutoff );
        }

        // Sort the lower structure
        std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
    }
    else if( nxSub >= nySub )
    {
        // Form the structure of a partition of the X dimension
        const int middle = (nxSub-1)/2;
        sn.size = nySub*nz;
        sn.offset = 
            ReorderedIndex
            ( xOffset+middle, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );

        // Allocate space for the lower structure
        int numJoins = 0;
        if( yOffset-1 >= 0 )
            ++numJoins;
        if( yOffset+nySub < ny )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        // Fill the (unsorted) lower structure
        int joinOffset = 0;
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset+middle, yOffset-1, i, nx, ny, nz, 
                  log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+middle, yOffset+nySub, i, nx, ny, nz, 
                  log2CommSize, cutoff );
        }

        // Sort the lower structure
        std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
    }
    else
    {
        // Form the structure of a partition of the Y dimension
        const int middle = (nySub-1)/2;
        sn.size = nxSub*nz;
        sn.offset = 
            ReorderedIndex
            ( xOffset, yOffset+middle, 0, nx, ny, nz, log2CommSize, cutoff );

        // Allocate space for the lower structure
        int numJoins = 0;
        if( xOffset-1 >= 0 )
            ++numJoins;
        if( xOffset+nxSub < nx )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        // Fill the (unsorted) lower structure
        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset-1, yOffset+middle, i, nx, ny, nz, 
                  log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+nxSub, yOffset+middle, i, nx, ny, nz,
                  log2CommSize, cutoff );
        }

        // Sort the lower structure
        std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
    }
#ifndef RELEASE
    clique::PopCallStack();
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

void FillLocalOrigStruct
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int log2CommSize, clique::symbolic::LocalSymmOrig& S )
{
#ifndef RELEASE
    clique::PushCallStack("FillLocalOrigStruct");
#endif
    // First count the depth, resize, and then run the actual fill
    int numSupernodes = 0;
    CountLocalTreeSize( nxSub, nySub, nz, cutoff, numSupernodes );
    S.supernodes.resize( numSupernodes );
    
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
    for( int s=numSupernodes-1; s>=0; --s )
    {
        Box box = boxStack.top();
        boxStack.pop();

        clique::symbolic::LocalSymmOrigSupernode& sn = S.supernodes[s];
        sn.parent = box.parentIndex;
        if( sn.parent != -1 )
        {
            if( box.leftChild )
                S.supernodes[sn.parent].children[0] = s;
            else
                S.supernodes[sn.parent].children[1] = s;
        }

        if( box.nx*box.ny*nz <= cutoff )
        {
            sn.size = box.nx*box.ny*nz;
            sn.offset = ReorderedIndex
                ( box.xOffset, box.yOffset, 0, nx, ny, nz, 
                  log2CommSize, cutoff );
            sn.children.clear();

            // Count, allocate, and fill the lower struct
            int numJoins = 0;
            if( box.xOffset-1 >= 0 )
                ++numJoins;
            if( box.xOffset+box.nx < nx )
                ++numJoins;
            if( box.yOffset-1 >= 0 )
                ++numJoins;
            if( box.yOffset+box.ny < ny )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            int joinOffset = 0;
            if( box.xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( box.xOffset-1, box.yOffset, i, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( box.xOffset+box.nx < nx )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( box.xOffset+box.nx, box.yOffset, i, nx, ny, nz,
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( box.yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( box.xOffset, box.yOffset-1, i, nx, ny, nz,
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( box.yOffset+box.ny < ny )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( box.xOffset, box.yOffset+box.ny, i, nx, ny, nz,
                      log2CommSize, cutoff );
            }

            // Sort the lower structure
            std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );
        }
        else
        {
            sn.children.resize(2);
            if( box.nx >= box.ny )
            {
                // Partition the X dimension (this is the separator)
                const int middle = (box.nx-1)/2;
                sn.size = box.ny*nz;
                sn.offset = ReorderedIndex
                    ( box.xOffset+middle, box.yOffset, 0, nx, ny, nz,
                      log2CommSize, cutoff );

                // Count, allocate, and fill the lower struct
                int numJoins = 0;
                if( box.yOffset-1 >= 0 )
                    ++numJoins;
                if( box.yOffset+box.ny < ny )
                    ++numJoins;
                sn.lowerStruct.resize( numJoins*nz );

                int joinOffset = 0;
                if( box.yOffset-1 >= 0 )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[i] = ReorderedIndex
                        ( box.xOffset+middle, box.yOffset-1, 0, nx, ny, nz,
                          log2CommSize, cutoff );
                    joinOffset += nz;
                }
                if( box.yOffset+box.ny < ny )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[joinOffset+i] = ReorderedIndex
                        ( box.xOffset+middle, box.yOffset+box.ny, 0, nx, ny, nz,
                          log2CommSize, cutoff );
                }

                // Sort the lower structure
                std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );

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
                sn.size = box.nx*nz;
                sn.offset = ReorderedIndex
                    ( box.xOffset, box.yOffset+middle, 0, nx, ny, nz,
                      log2CommSize, cutoff );

                // Count, allocate, and fill the lower struct
                int numJoins = 0;
                if( box.xOffset-1 >= 0 )
                    ++numJoins;
                if( box.xOffset+box.nx < nx )
                    ++numJoins;
                sn.lowerStruct.resize( numJoins*nz );

                int joinOffset = 0;
                if( box.xOffset-1 >= 0 )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[i] = ReorderedIndex
                        ( box.xOffset-1, box.yOffset+middle, 0, nx, ny, nz,
                          log2CommSize, cutoff );
                    joinOffset += nz;
                }
                if( box.xOffset+box.nx < nx )
                {
                    for( int i=0; i<nz; ++i )
                        sn.lowerStruct[joinOffset+i] = ReorderedIndex
                        ( box.xOffset+box.nx, box.yOffset+middle, 0, nx, ny, nz,
                          log2CommSize, cutoff );
                }

                // Sort the lower structure
                std::sort( sn.lowerStruct.begin(), sn.lowerStruct.end() );

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
    clique::PopCallStack();
#endif
}

void CountLocalTreeSize
( int nx, int ny, int nz, int cutoff, int& numSupernodes )
{
    ++numSupernodes;
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
        CountLocalTreeSize( middle, ny, nz, cutoff, numSupernodes );

        // Recurse on the right partition
        CountLocalTreeSize
        ( std::max(nx-middle-1,0), ny, nz, cutoff, numSupernodes );
    }
    else
    {
        // Partition the Y dimension
        const int middle = (ny-1)/2;

        // Recurse on the top partition
        CountLocalTreeSize( nx, middle, nz, cutoff, numSupernodes );

        // Recurse on the bottom partition
        CountLocalTreeSize
        ( nx, std::max(ny-middle-1,0), nz, cutoff, numSupernodes );
    }
}

