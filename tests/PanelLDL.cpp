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

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

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

    // Fill the distributed portion of the structure
    int nxSub=nx, nySub=ny, xOffset=0, yOffset=0;
    clique::symbolic::DistSymmOrig distOrig;
    distOrig.supernodes.resize( log2CommSize );
    for( int s=0; s<log2CommSize; ++s )
    {
        clique::symbolic::DistSymmOrigSupernode& sn = 
            distOrig.supernodes[log2CommSize-(s+1)];
        if( nxSub >= nySub )
        {
            // Form the structure of a partition of the X dimension
            sn.size = nySub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );

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
                    ( xOffset, yOffset-1, 0, nx, ny, nz, log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( yOffset+nySub < ny )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset, yOffset+nySub, 0, nx, ny, nz, 
                      log2CommSize, cutoff );
            }
            
            // Pick the new offsets and sizes based upon whether or not the
            // (log2CommSize-s)'th bit of our rank is nonzero.
            // HERE
        }
        else
        {
            // Form the structure of a partition of the Y dimension
            sn.size = nxSub*nz;
            sn.offset = 
                ReorderedIndex
                ( xOffset, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );

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
                    ( xOffset-1, yOffset, 0, nx, ny, nz, log2CommSize, cutoff );
                joinOffset += nz;
            }
            if( xOffset+nxSub < nx )
            {
                for( int i=0; i<nz; ++i )
                    sn.lowerStruct[joinOffset+i] = ReorderedIndex
                    ( xOffset+nxSub, yOffset, 0, nx, ny, nz,
                      log2CommSize, cutoff );
            }

            // Pick the new offsets and sizes based upon whether or not the
            // (log2CommSize-s)'th bit of our rank is nonzero.
            // HERE
        }
    }
    // TODO: Form the local portion of the structure.

    // TODO: Call the symbolic factorization routine
    // TODO: Call the numerical factorization routine
    // TODO: Call a solve routine

    clique::Finalize();
    return 0;
}

