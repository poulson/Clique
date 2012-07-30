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
      << "NestedDissection <n> [cutoff=128] [numDistSeps=10] [numSeqSeps=5]\n" 
      << "  n: size of n x n x n mesh\n"
      << "  cutoff: maximum size of leaf node\n"
      << "  numDistSeps: number of distributed separators to try\n"
      << "  numSeqSeps: number of sequential separators to try\n"
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
    const int cutoff = ( argc >= 3 ? atoi( argv[2] ) : 128 );
    const int numDistSeps = ( argc >= 4 ? atoi( argv[3] ) : 10 );
    const int numSeqSeps = ( argc >= 5 ? atoi( argv[4] ) : 5 );

    try
    {
        const int numVertices = n*n*n;
        DistGraph graph( numVertices, comm );

        const int firstLocalSource = graph.FirstLocalSource();
        const int numLocalSources = graph.NumLocalSources();

        // Fill our portion of the graph of a 3D n x n x n 7-point stencil
        // in natural ordering: (x,y,z) at x + y*n + z*n*n
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of graph...";
            std::cout.flush();
        }
        graph.StartAssembly();
        graph.Reserve( 7*numLocalSources );
        for( int iLocal=0; iLocal<numLocalSources; ++iLocal )
        {
            const int i = firstLocalSource + iLocal;
            const int x = i % n;
            const int y = (i/n) % n;
            const int z = i/(n*n);

            graph.Insert( i, i );
            if( x != 0 )
                graph.Insert( i, i-1 );
            if( x != n-1 )
                graph.Insert( i, i+1 );
            if( y != 0 )
                graph.Insert( i, i-n );
            if( y != n-1 )
                graph.Insert( i, i+n );
            if( z != 0 )
                graph.Insert( i, i-n*n );
            if( z != n-1 )
                graph.Insert( i, i+n*n );
        }
        graph.StopAssembly();
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running nested dissection...";
            std::cout.flush();
        }
        DistSymmInfo info;
        DistSeparatorTree sepTree;
        DistMap map;
        NestedDissection
        ( graph, map, sepTree, info, cutoff, numDistSeps, numSeqSeps );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "done" << std::endl;

        if( commRank == 0 )
        {
            const int numDistNodes = info.distNodes.size();
            const int numLocalNodes = info.localNodes.size(); 
            const int rootSepSize = info.distNodes.back().size;
            std::cout << "\n"
                      << "On the root process:\n"
                      << "-----------------------------------------\n"
                      << numLocalNodes << " local nodes\n"
                      << numDistNodes  << " distributed nodes\n"
                      << rootSepSize << " vertices in root separator\n"
                      << "\n";
            for( int s=0; s<rootSepSize; ++s )
            {
                const int i = 
                    ( numDistNodes > 1 ? 
                      sepTree.distSeps.back().indices[s] :
                      sepTree.localSepsAndLeaves.back()->indices[s] );
                const int x = i % n;
                const int y = (i/n) % n;
                const int z = i/(n*n);
                std::cout << "rootSep[" << s << "]: " << i << ", ("
                          << x << "," << y << "," << z << ")\n";
            }
            std::cout << std::endl;
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
