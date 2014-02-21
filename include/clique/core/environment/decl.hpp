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

using elem::byte;
using elem::Int;
using elem::Unsigned;

// Pull in some of Elemental's imported libraries
namespace blas = elem::blas;
namespace lapack = elem::lapack;
namespace mpi = elem::mpi;

// Pull in a number of useful enums from Elemental
using namespace elem::distribution_wrapper;
using namespace elem::left_or_right_wrapper;
using namespace elem::orientation_wrapper;
using namespace elem::unit_or_non_unit_wrapper;
using namespace elem::upper_or_lower_wrapper;

// For scalar operations
using elem::Base;
using elem::Complex;
using elem::Abs;
using elem::Sqrt;

// Pull in a few classes from Elemental
using elem::Matrix;
using elem::Grid;
using elem::DistMatrix;
using elem::View;
using elem::LockedView;

// Pull in a few indexing routines
using elem::Shift;
using elem::Length;

// For easily throwing errors
using elem::LogicError;
using elem::RuntimeError;

using elem::SwapClear;

void PrintVersion( std::ostream& os=std::cout );
void PrintConfig( std::ostream& os=std::cout );
void PrintCCompilerInfo( std::ostream& os=std::cout );
void PrintCxxCompilerInfo( std::ostream& os=std::cout );
 
bool Initialized();
void Initialize( int& argc, char**& argv );
void Finalize();

// For getting the MPI argument instance (for internal usage)
class Args : public elem::choice::MpiArgs
{
public:
    Args
    ( int argc, char** argv,
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr )
    : elem::choice::MpiArgs(argc,argv,comm,error)
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
