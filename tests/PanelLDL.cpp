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
    std::cout << "PanelLDL <nx> <ny> <nz>\n"
              << "<nx>: size of panel in x direction\n"
              << "<ny>: size of panel in y direction\n"
              << "<nz>: size of panel in z direction\n" << std::endl;
}

int ReorderedIndex
( int nx, int ny, int nz, int x, int y, int z, 
  int minBalancedDepth, int maxLeafSize )
{
    // HERE: Call a different recursive routine?
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

    if( argc < 4 )
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
    if( commRank == 0 )
        std::cout << "(nx,ny,nz)=(" << nx << "," << ny << "," << nz << ")"
                  << std::endl;

    // TODO: Perform the analytical nested dissection and fill the distributed
    //       then local original structures.
    clique::symbolic::DistSymmOrig distSymmOrig;
    clique::symbolic::LocalSymmOrig localSymmOrig;
    for( int s=0; s<log2CommSize; ++s )
    {
        // HERE
    }
    // TODO: Call the symbolic factorization routine
    // TODO: Call the numerical factorization routine
    // TODO: Call a solve routine

    clique::Finalize();
    return 0;
}

