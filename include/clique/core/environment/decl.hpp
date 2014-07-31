/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_CORE_ENVIRONMENT_DECL_HPP
#define CLIQ_CORE_ENVIRONMENT_DECL_HPP

namespace cliq {

using El::byte;
using El::Int;
using El::Unsigned;

// Pull in some of Elemental's imported libraries
namespace blas = El::blas;
namespace lapack = El::lapack;
namespace mpi = El::mpi;

// Pull in a number of useful enums from Elemental
using namespace El::DistNS;
using namespace El::LeftOrRightNS;
using namespace El::OrientationNS;
using namespace El::UnitOrNonUnitNS;
using namespace El::UpperOrLowerNS;

// For scalar operations
using El::Base;
using El::Complex;
using El::Abs;
using El::Sqrt;

// Pull in a few classes from Elemental
using El::Matrix;
using El::Grid;
using El::DistMatrix;
using El::View;
using El::LockedView;

// Pull in a few indexing routines
using El::Shift;
using El::Length;

// For easily throwing errors
using El::LogicError;
using El::RuntimeError;

using El::SwapClear;

void PrintVersion( std::ostream& os=std::cout );
void PrintConfig( std::ostream& os=std::cout );
void PrintCCompilerInfo( std::ostream& os=std::cout );
void PrintCxxCompilerInfo( std::ostream& os=std::cout );
 
bool Initialized();
void Initialize( int& argc, char**& argv );
void Finalize();

// For getting the MPI argument instance (for internal usage)
class Args : public El::choice::MpiArgs
{
public:
    Args
    ( int argc, char** argv,
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr )
    : El::choice::MpiArgs(argc,argv,comm,error)
    { }
    virtual ~Args() { }
protected:
    virtual void HandleVersion( std::ostream& os=std::cout ) const;
    virtual void HandleBuild( std::ostream& os=std::cout ) const;
};
Args& GetArgs();

// For processing command-line arguments
template<typename T>
T Input( std::string name, std::string desc );
template<typename T>
T Input( std::string name, std::string desc, T defaultVal );
void ProcessInput();
void PrintInputReport();

DEBUG_ONLY(
    void PushCallStack( std::string s );
    void PopCallStack();
    void DumpCallStack();

    class CallStackEntry
    {
    public:
        CallStackEntry( std::string s ) 
        { 
            if( !std::uncaught_exception() )
                PushCallStack(s);
        }
        ~CallStackEntry() 
        {     
            if( !std::uncaught_exception() )
                PopCallStack(); 
        }
    };
)

void ReportException( std::exception& e );

template<typename T>
struct Entry
{
    Int i, j;
    T value;
};

template<typename T>
bool IsSorted( const std::vector<T>& x );
// While is_strictly_sorted exists in Boost, it does not exist in the STL (yet)
template<typename T>
bool IsStrictlySorted( const std::vector<T>& x );

void Union
( std::vector<Int>& both, 
  const std::vector<Int>& first, const std::vector<Int>& second );
std::vector<Int> 
Union( const std::vector<Int>& first, const std::vector<Int>& second );

void RelativeIndices
( std::vector<Int>& relInds, 
  const std::vector<Int>& sub, const std::vector<Int>& full );
std::vector<Int>
RelativeIndices( const std::vector<Int>& sub, const std::vector<Int>& full );

Int RowToProcess( Int i, Int blocksize, Int commSize );

Int Find
( const std::vector<Int>& sortedInds, Int index, 
  std::string msg="Could not find index" );

void
VerifySendsAndRecvs
( const std::vector<int>& sendCounts,
  const std::vector<int>& recvCounts, mpi::Comm comm );

template<typename T>
void SparseAllToAll
( const std::vector<T>& sendBuffer,
  const std::vector<int>& sendCounts, const std::vector<int>& sendDispls,
        std::vector<T>& recvBuffer,
  const std::vector<int>& recvCounts, const std::vector<int>& recvDispls,
        mpi::Comm comm );

} // namespace cliq

#endif // ifndef CLIQ_CORE_ENVIRONMENT_DECL_HPP
