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

namespace { bool cliqueInitializedMpi; }

void clique::Initialize( int& argc, char**& argv )
{
    int mpiInitialized;
    MPI_Initialized( &mpiInitialized );
    if( !mpiInitialized )
    {
        int mpiFinalized;    
        MPI_Finalized( &mpiFinalized );
        if( mpiFinalized )
            throw std::logic_error
            ("Cannot initialize Clique after finalizing MPI");
#ifndef _OPENMP
        int provided;
        MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
        if( provided != MPI_THREAD_MULTIPLE )
            std::cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                      << std::endl;
#else
        MPI_Init( &argc, &argv );
#endif // ifndef _OPENMP
        ::cliqueInitializedMpi = true;
    }
    else
    {
        ::cliqueInitializedMpi = false;
    }
}

void clique::Finalize()
{
    int mpiFinalized;
    MPI_Finalized( &mpiFinalized );
    if( mpiFinalized )
        std::cerr << "WARNING: MPI was finalized before Clique." << std::endl;
    else if( ::cliqueInitializedMpi )
        MPI_Finalize();
}

#ifndef RELEASE
namespace { std::stack<std::string> callStack; }

void clique::PushCallStack( std::string s )
{ ::callStack.push( s ); }

void clique::PopCallStack()
{ ::callStack.pop(); }

void clique::DumpCallStack()
{
    std::ostringstream msg;
    while( !::callStack.empty() )
    {
        msg << "[" << ::callStack.size() << "]: " << ::callStack.top() << "\n";
        ::callStack.pop();
    }
    std::cerr << msg.str() << std::endl;
}
#endif // ifndef RELEASE
