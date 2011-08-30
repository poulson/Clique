/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
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

void Usage()
{
    std::cout << 
        "DistDenseSymmMatrix <gridHeight> <gridWidth> <n> <blockSize>\n" 
        << std::endl;
}

int
main( int argc, char* argv[] )
{
    clique::Initialize( argc, argv );
    clique::mpi::Comm comm = clique::mpi::COMM_WORLD;
    const int commRank = clique::mpi::CommRank( comm );
    const int commSize = clique::mpi::CommSize( comm );

    if( argc < 5 )
    {
        if( commRank == 0 ) 
            Usage();
        clique::Finalize();
        return 0;
    }

    const int gridHeight = atoi(argv[1]);
    const int gridWidth = atoi(argv[2]);
    const int height = atoi(argv[3]);
    const int blockSize = atoi(argv[4]);

    if( gridHeight*gridWidth != commSize )
    {
        if( commRank == 0 )
            std::cerr << "Grid dimensions do not match number of processes"
                      << std::endl;
        clique::Finalize();
        return 0;
    }

    clique::DistDenseSymmMatrix<double> A( comm, gridHeight, gridWidth );
    A.Print("Empty A");
    A.Reconfigure( height, blockSize );
    A.Print("Default reconfigured A");
    A.MakeIdentity();
    A.Print("A := I");

    clique::Finalize();
    return 0;
}

