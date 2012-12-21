/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
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
    std::cerr << msg.str() << std::endl;
}
#endif // ifndef RELEASE

} // namespace cliq
