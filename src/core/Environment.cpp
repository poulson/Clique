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
    if( !mpi::Initialized() )
    {
        if( mpi::Finalized() )
            throw std::logic_error
            ("Cannot initialize Clique after finalizing MPI");
#ifndef _OPENMP
        const int provided = 
            mpi::InitThread( argc, argv, mpi::THREAD_MULTIPLE );
        if( provided != mpi::THREAD_MULTIPLE )
            std::cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                      << std::endl;
#else
        mpi::Init( argc, argv );
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
    if( mpi::Finalized() )
        std::cerr << "WARNING: MPI was finalized before Clique." << std::endl;
    else if( ::cliqueInitializedMpi )
        mpi::Finalize();
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
