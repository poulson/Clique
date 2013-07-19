/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "clique.hpp"

namespace { 
bool cliqueInitializedElemental; 
bool initializedClique = false;
#ifndef RELEASE
std::stack<std::string> callStack;
#endif
}

namespace cliq {

void PrintVersion( std::ostream& os )
{
    os << "Clique version information:\n"
       << "  Git revision: " << GIT_SHA1 << "\n"
       << "  Version:      " << Clique_VERSION_MAJOR << "."
                             << Clique_VERSION_MINOR << "\n"
       << "  Build type:   " << CMAKE_BUILD_TYPE << "\n"
       << std::endl;
}

void PrintCCompilerInfo( std::ostream& os )
{
    os << "Clique's C compiler info:\n"
       << "  CMAKE_C_COMPILER:    " << CMAKE_C_COMPILER << "\n"
       << "  MPI_C_COMPILER:      " << MPI_C_COMPILER << "\n"
       << "  MPI_C_INCLUDE_PATH:  " << MPI_C_INCLUDE_PATH << "\n"
       << "  MPI_C_COMPILE_FLAGS: " << MPI_C_COMPILE_FLAGS << "\n"
       << "  MPI_C_LINK_FLAGS:    " << MPI_C_LINK_FLAGS << "\n"
       << "  MPI_C_LIBRARIES:     " << MPI_C_LIBRARIES << "\n"
       << std::endl;
}

void PrintCxxCompilerInfo( std::ostream& os )
{
    os << "Clique's C++ compiler info:\n"
       << "  CMAKE_CXX_COMPILER:    " << CMAKE_CXX_COMPILER << "\n"
       << "  CXX_FLAGS:             " << CXX_FLAGS << "\n"
       << "  MPI_CXX_COMPILER:      " << MPI_CXX_COMPILER << "\n"
       << "  MPI_CXX_INCLUDE_PATH:  " << MPI_CXX_INCLUDE_PATH << "\n"
       << "  MPI_CXX_COMPILE_FLAGS: " << MPI_CXX_COMPILE_FLAGS << "\n"
       << "  MPI_CXX_LINK_FLAGS:    " << MPI_CXX_LINK_FLAGS << "\n"
       << "  MPI_CXX_LIBRARIES:     " << MPI_CXX_LIBRARIES << "\n"
       << std::endl;
}

bool Initialized()
{ return ::initializedClique; }

void Initialize( int& argc, char**& argv )
{
    // If Clique has already been initialized, this is a no-op
    if( ::initializedClique )
        return;

    const bool mustInitElemental = !elem::Initialized();
    if( mustInitElemental )
    {
        elem::Initialize( argc, argv );
        ::cliqueInitializedElemental = true;
    }
    else
    {
        ::cliqueInitializedElemental = false;
    }
    ::initializedClique = true;
}

void Finalize()
{
    // If Clique is not currently initialized, then this is a no-op
    if( !::initializedClique )
        return;
    
    if( ::cliqueInitializedElemental )
        elem::Finalize();

    ::initializedClique = false;
}

#ifndef RELEASE
void PushCallStack( std::string s )
{ ::callStack.push( s ); }

void PopCallStack()
{ ::callStack.pop(); }

void DumpCallStack()
{
    std::ostringstream msg;
    while( !::callStack.empty() )
    {
        msg << "[" << ::callStack.size() << "]: " << ::callStack.top() << "\n";
        ::callStack.pop();
    }
    std::cerr << msg.str();;
    std::cerr.flush();
}
#endif // ifndef RELEASE

void ReportException( std::exception& e )
{
    elem::ReportException( e );
#ifndef RELEASE
    DumpCallStack();
#endif
}

} // namespace cliq
