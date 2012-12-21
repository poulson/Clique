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

    try
    {
        const unsigned numMB = Input
            ("--numMB","number of megabytes to read/write",100u);
        const unsigned numPieces = Input
            ("--numPieces","number of pieces to break transfer into",5);
        const unsigned numFiles = Input
            ("--numFiles","number of files to read/write",5);
        const std::string baseName = Input
            ("--baseName","base name for files",std::string("scratch"));
        ProcessInput();

        if( numPieces == 0 && commRank == 0 )
        {
            std::cout << "Number of pieces must be positive." << std::endl;
            cliq::Finalize();
            return 0;
        }
        if( numFiles == 0 && commRank == 0 )
        {
            std::cout << "Number of files must be positive." << std::endl;
            cliq::Finalize();
            return 0;
        }
        if( commRank == 0 )
            std::cout << "numMB:     " << numMB << "\n"
                      << "numPieces: " << numPieces << "\n"
                      << "numFiles:  " << numFiles << "\n"
                      << "baseName:  " << baseName << "\n"
                      << std::endl;

        const std::size_t bufferSize = numMB<<20;
        std::vector<char> buffer( bufferSize );

        const std::size_t normalPieceSize = bufferSize / numPieces;
        const std::size_t lastPieceSize = 
            normalPieceSize + (bufferSize%numPieces);

        // Write the files
        for( unsigned f=0; f<numFiles; ++f )
        {
            std::ostringstream filename;
            filename << baseName << "-" << commRank << "-" << f << ".dat";

            mpi::Barrier( comm );
            const double writeStartTime = mpi::Time();
            std::ofstream outFile
            ( filename.str().c_str(), std::ios::out|std::ios::binary );
            for( unsigned i=0; i<numPieces-1; ++i )
                outFile.write( &buffer[i*normalPieceSize], normalPieceSize );
            outFile.write
            ( &buffer[(numPieces-1)*normalPieceSize], lastPieceSize );
            outFile.close();
            mpi::Barrier( comm );
            const double writeStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Write time: " << writeStopTime-writeStartTime
                          << " secs." << std::endl;
        }

        // Read the files
        for( unsigned f=0; f<numFiles; ++f )
        {
            std::ostringstream filename;
            filename << baseName << "-" << commRank << "-" << f << ".dat";

            mpi::Barrier( comm );
            const double readStartTime = mpi::Time();
            std::ifstream inFile
            ( filename.str().c_str(), std::ios::in|std::ios::binary );
            for( unsigned i=0; i<numPieces-1; ++i )
            inFile.read( &buffer[i*normalPieceSize], normalPieceSize );
            inFile.read
            ( &buffer[(numPieces-1)*normalPieceSize], lastPieceSize );
            inFile.close();
            mpi::Barrier( comm );
            const double readStopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << "Read time: " << readStopTime-readStartTime
                          << " secs." << std::endl;
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
