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
( int nxSub, int nySub, int nz, int xOffset, int yOffset, int cutoff, 
  int& numSupernodes );

void
FillOrigStruct
( int nx, int ny, int nz, int cutoff, int commRank, int log2CommSize,
  clique::symbolic::LocalSymmOrig& localOrig, 
  clique::symbolic::DistSymmOrig& distOrig );

void
FillDistOrigStruct
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, int commRank, int log2CommSize,
  clique::symbolic::DistSymmOrig& SOrig );

void
FillLocalOrigStruct
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

    // Fill the distributed portion of the original structure
    clique::symbolic::DistSymmOrig distSOrig;
    clique::symbolic::LocalSymmOrig localSOrig;
    FillOrigStruct
    ( nx, ny, nz, cutoff, commRank, log2CommSize, localSOrig, distSOrig );

    // Call the symbolic factorization routine
    clique::symbolic::DistSymmFact distS;
    clique::symbolic::LocalSymmFact localS;
    clique::symbolic::SymmetricFactorization
    ( localSOrig, distSOrig, localS, distS, true );

    // Directly initialize the frontal matrices with the original sparse matrix
    clique::numeric::DistSymmFact<F> distL;
    clique::numeric::LocalSymmFact<F> localL;
    // TODO

    // Call the numerical factorization routine
    //clique::numeric::LDL( ADJOINT, localS, distS, localL, distL );

    // Set up the properly ordered RHS and call a solve routine
    //clique::numeric::LDLSolve
    //( ADJOINT, localS, distS, localL, distL, localX, true );

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
    PushCallStack("ReorderedIndex");
#endif
    int index = ReorderedIndexRecursion
    ( x, y, z, nx, ny, nz, minBalancedDepth, cutoff, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
    return index;
}

void
FillOrigStruct
( int nx, int ny, int nz, int cutoff, int commRank, int log2CommSize,
  clique::symbolic::LocalSymmOrig& localSOrig, 
  clique::symbolic::DistSymmOrig& distSOrig )
{
    int nxSub=nx, nySub=ny, xOffset=0, yOffset=0;
    FillDistOrigStruct
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, 
      cutoff, commRank, log2CommSize, distSOrig );
    FillLocalOrigStruct
    ( nx, ny, nz, nxSub, nySub, xOffset, yOffset, cutoff, 
      log2CommSize, localSOrig );
}
  
void
FillDistOrigStruct
( int nx, int ny, int nz, int& nxSub, int& nySub, int& xOffset, int& yOffset, 
  int cutoff, int commRank, int log2CommSize, 
  clique::symbolic::DistSymmOrig& SOrig )
{
    SOrig.supernodes.resize( log2CommSize );
    // Fill the distributed nodes
    for( int s=log2CommSize-1; s>0; --s )
    {
        clique::symbolic::DistSymmOrigSupernode& sn = SOrig.supernodes[s];
        const bool powerOfTwo = 1u<<(s-1);
        const bool onLeft = (commRank&powerOfTwo) == 0;
        if( nxSub >= nySub )
        {
            // Form the structure of a partition of the X dimension
            sn.size = nySub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );
            const int middle = (nxSub-1)/2;

            int numJoins = 0;
            if( yOffset-1 >= 0 )
                ++numJoins;
            if( yOffset+nySub < ny )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            int joinOffset = 0;
            if( yOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( xOffset+middle, yOffset-1, 0, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( yOffset+nySub < ny )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+middle, yOffset+nySub, 0, nx, ny, nz, 
                      log2CommSize, cutoff );
            }
            
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
            sn.size = nxSub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );
            const int middle = (nySub-1)/2;

            int numJoins = 0;
            if( xOffset-1 >= 0 )
                ++numJoins;
            if( xOffset+nxSub < nx )
                ++numJoins;
            sn.lowerStruct.resize( numJoins*nz );

            int joinOffset = 0;
            if( xOffset-1 >= 0 )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[i] = ReorderedIndex
                    ( xOffset-1, yOffset+middle, 0, nx, ny, nz, 
                      log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( xOffset+nxSub < nx )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+nxSub, yOffset+middle, 0, nx, ny, nz,
                      log2CommSize, cutoff );
            }

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
    if( nxSub >= nySub )
    {
        // Form the structure of a partition of the X dimension
        sn.size = nySub*nz;
        sn.offset = 
            ReorderedIndex
            ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );
        const int middle = (nxSub-1)/2;

        int numJoins = 0;
        if( yOffset-1 >= 0 )
            ++numJoins;
        if( yOffset+nySub < ny )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        int joinOffset = 0;
        if( yOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset+middle, yOffset-1, 0, nx, ny, nz, 
                  log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( yOffset+nySub < ny )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+middle, yOffset+nySub, 0, nx, ny, nz, 
                  log2CommSize, cutoff );
        }
    }
    else
    {
        // Form the structure of a partition of the Y dimension
        sn.size = nxSub*nz;
        sn.offset = 
            ReorderedIndex
            ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );
        const int middle = (nySub-1)/2;

        int numJoins = 0;
        if( xOffset-1 >= 0 )
            ++numJoins;
        if( xOffset+nxSub < nx )
            ++numJoins;
        sn.lowerStruct.resize( numJoins*nz );

        int joinOffset = 0;
        if( xOffset-1 >= 0 )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[i] = ReorderedIndex
                ( xOffset-1, yOffset+middle, 0, nx, ny, nz, 
                  log2CommSize, cutoff );
            joinOffset += nz;
        }
        if( xOffset+nxSub < nx )
        {
            for( int i=0; i<nz; ++i )
                sn.lowerStruct[joinOffset+i] = ReorderedIndex
                ( xOffset+nxSub, yOffset+middle, 0, nx, ny, nz,
                  log2CommSize, cutoff );
        }
    }
}

// Even though this is only used within the following function, it must be
// declared outside of the function in order to be used as a template argument
// (for std::stack) since template arguments must have external linkage.
//
// This data holds the bounding box information.
struct Box
{
    int parentIndex, nx, ny, xOffset, yOffset;
    bool marked;
};

void
FillLocalOrigStruct
( int nx, int ny, int nz, int nxSub, int nySub, int xOffset, int yOffset, 
  int cutoff, int log2CommSize, clique::symbolic::LocalSymmOrig& S )
{
    // First count the depth, resize, and then run the actual fill
    int numSupernodes = 0;
    CountLocalTreeSize
    ( nxSub, nySub, nz, xOffset, yOffset, cutoff, numSupernodes );
    S.supernodes.resize( numSupernodes );
    std::stack<Box> boxStack;

    // Initialize the root's box
    {
        Box box;
        box.parentIndex = -1;
        box.nx = nxSub;
        box.ny = nySub;
        box.xOffset = xOffset;
        box.yOffset = yOffset;
        box.marked = false;
        boxStack.push(box);
    }

    for( int s=numSupernodes-1; s>=0; --s )
    {
        Box box = boxStack.top();

        clique::symbolic::LocalSymmOrigSupernode& sn = S.supernodes[s];
        sn.parent = box.parentIndex;
        if( sn.parent != -1 )
        {
            if( box.marked )
                S.supernodes[sn.parent].children[1] = s;
            else
                S.supernodes[sn.parent].children[0] = s;
        }

        if( box.nx*box.ny*nz <= cutoff )
        {
            sn.size = box.nx*box.ny*nz;
            sn.offset = ReorderedIndex
                ( box.xOffset, box.yOffset, 0, nx, ny, nz, 
                  log2CommSize, cutoff );
            sn.children.clear();

            // Count, allocate, and fill the lower struct
            // Essentially, we need to check all four quasi2d borders
            // TODO
            
            boxStack.pop();
            if( !boxStack.empty() )
                boxStack.top().marked = true;
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

                if( box.marked )
                {
                    // Set up for the right child
                    Box rightBox;
                    rightBox.parentIndex = s;
                    rightBox.nx = std::max(box.nx-middle-1,0);
                    rightBox.ny = box.ny;
                    rightBox.xOffset = box.xOffset+middle+1;
                    rightBox.yOffset = box.yOffset;
                    rightBox.marked = false;
                    boxStack.push( rightBox );
                }
                else
                {
                    // Set up for the left child
                    Box leftBox;
                    leftBox.parentIndex = s;
                    leftBox.nx = middle;
                    leftBox.ny = box.ny;
                    leftBox.xOffset = box.xOffset;
                    leftBox.yOffset = box.yOffset;
                    leftBox.marked = false;
                    boxStack.push( leftBox );
                }
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

                if( box.marked )
                {
                    // Set up for the right child
                    Box rightBox;
                    rightBox.parentIndex = s;
                    rightBox.nx = box.nx;
                    rightBox.ny = std::max(box.ny-middle-1,0);
                    rightBox.xOffset = box.xOffset;
                    rightBox.yOffset = box.yOffset+middle+1;
                    rightBox.marked = false;
                    boxStack.push( rightBox );
                }
                else
                {
                    // Set up for the left child
                    Box leftBox;
                    leftBox.parentIndex = s;
                    leftBox.nx = box.nx;
                    leftBox.ny = middle;
                    leftBox.xOffset = box.xOffset;
                    leftBox.yOffset = box.yOffset;
                    leftBox.marked = false;
                    boxStack.push( leftBox );
                }
            }
        }
    }
}

void CountLocalTreeSize
( int nx, int ny, int nz, int xOffset, int yOffset, int cutoff, 
  int& numSupernodes )
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
        CountLocalTreeSize
        ( middle, ny, nz, xOffset, yOffset, cutoff, numSupernodes );

        // Recurse on the right partition
        CountLocalTreeSize
        ( std::max(nx-middle-1,0), ny, nz, xOffset+middle+1, yOffset, 
          cutoff, numSupernodes );
    }
    else
    {
        // Partition the Y dimension
        const int middle = (ny-1)/2;

        // Recurse on the top partition
        CountLocalTreeSize
        ( nx, middle, nz, xOffset, yOffset, cutoff, numSupernodes );

        // Recurse on the bottom partition
        CountLocalTreeSize
        ( nx, std::max(ny-middle-1,0), nz, xOffset, yOffset+middle+1, 
          cutoff, numSupernodes );
    }
}

