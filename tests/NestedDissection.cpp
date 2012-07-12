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
    std::cout << "NestedDissection <n>\n" 
              << "  n: size of n x n x n mesh\n"
              << std::endl;
}

int
main( int argc, char* argv[] )
{
    cliq::Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 2 )
    {
        if( commRank == 0 )
            Usage();
        cliq::Finalize();
        return 0;
    }
    const int n = atoi( argv[1] );

    try
    {
        const int numVertices = n*n*n;
        DistGraph graph( numVertices, comm );

        const int firstLocalSource = graph.FirstLocalSource();
        const int numLocalSources = graph.NumLocalSources();

        // Fill our portion of the graph of a 3D n x n x n 7-point stencil
        // in natural ordering: (x,y,z) at x + y*n + z*n*n
        graph.StartAssembly();
        graph.Reserve( 7*numLocalSources );
        for( int iLocal=0; iLocal<numLocalSources; ++iLocal )
        {
            const int i = firstLocalSource + iLocal;
            const int x = i % n;
            const int y = (i/n) % n;
            const int z = i/(n*n);

            if( z != 0 )
                graph.PushBack( i, i-n*n );
            if( y != 0 )
                graph.PushBack( i, i-n );
            if( x != 0 )
                graph.PushBack( i, i-1 );
            graph.PushBack( i, i );
            if( x != n-1 )
                graph.PushBack( i, i+1 );
            if( y != n-1 )
                graph.PushBack( i, i+n );
            if( z != n-1 )
                graph.PushBack( i, i+n*n );
        }
        graph.StopAssembly();

        DistGraph child;
        std::vector<int> localMap;
        bool haveLeftChild;
        const int sepSize = Bisect( graph, child, localMap, haveLeftChild );

        int leftSize, rightSize;
        if( haveLeftChild )
        {
            leftSize = child.NumSources();
            rightSize = numVertices - leftSize - sepSize;
        }
        else
        {
            rightSize = child.NumSources();
            leftSize = numVertices - rightSize - sepSize;
        }
        if( commRank == 0 )
        {
            if( haveLeftChild )
                std::cout << "Root is on left with sizes: " 
                          << leftSize << ", " << rightSize << ", " << sepSize
                          << std::endl;
            else
                std::cout << "Root is on right with sizes: " 
                          << leftSize << ", " << rightSize << ", " << sepSize
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
