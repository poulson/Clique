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

namespace {

inline void 
SafeMpi( int mpiError )
{
#ifndef RELEASE
    if( mpiError != MPI_SUCCESS )    
    {
        char errorString[200];
        int lengthOfErrorString;
        MPI_Error_string( mpiError, errorString, &lengthOfErrorString );
        throw std::logic_error( errorString );
    }
#endif
}

} // anonymous namespace

//----------------------------//
// MPI environmental routines //
//----------------------------//

void
clique::mpi::Init( int& argc, char**& argv )
{ MPI_Init( &argc, &argv ); }

int 
clique::mpi::InitThread( int& argc, char**& argv, int required )
{ 
    int provided; 
    MPI_Init_thread( &argc, &argv, required, &provided ); 
    return provided;
}

void
clique::mpi::Finalize()
{ MPI_Finalize(); }

int
clique::mpi::Initialized()
{ 
    int initialized;
    MPI_Initialized( &initialized );
    return initialized;
}

int
clique::mpi::Finalized()
{
    int finalized;
    MPI_Finalized( &finalized );
    return finalized;
}

double
clique::mpi::Time()
{ return MPI_Wtime(); }


void
clique::mpi::OpCreate
( mpi::UserFunction* func, int commutes, mpi::Op& op )
{
#ifndef RELEASE
    PushCallStack("mpi::OpCreate");
#endif
    SafeMpi( MPI_Op_create( func, commutes, &op ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::OpFree( mpi::Op& op )
{
#ifndef RELEASE
    PushCallStack("mpi::OpFree");
#endif
    SafeMpi( MPI_Op_free( &op ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//---------------------------//
// Communicator manipulation //
//---------------------------//

int
clique::mpi::CommRank( mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommRank");
#endif
    int rank;
    SafeMpi( MPI_Comm_rank( comm, &rank ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

int
clique::mpi::CommSize( mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommSize");
#endif
    int size;
    SafeMpi( MPI_Comm_size( comm, &size ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return size;
}

void
clique::mpi::CommCreate
( mpi::Comm parentComm, mpi::Group subsetGroup, mpi::Comm& subsetComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommCreate");
#endif
    SafeMpi( MPI_Comm_create( parentComm, subsetGroup, &subsetComm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::CommDup( mpi::Comm original, mpi::Comm& duplicate )
{
#ifndef RELEASE
    PushCallStack("mpi::CommDup");
#endif
    SafeMpi( MPI_Comm_dup( original, &duplicate ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::CommSplit
( mpi::Comm comm, int color, int key, mpi::Comm& newComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommSplit");
#endif
    SafeMpi( MPI_Comm_split( comm, color, key, &newComm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::CommFree( mpi::Comm& comm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommFree");
#endif
    SafeMpi( MPI_Comm_free( &comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

bool
clique::mpi::CongruentComms( mpi::Comm comm1, mpi::Comm comm2 )
{
#ifndef RELEASE
    PushCallStack("mpi::CongruentComms");
#endif
    int result;
    SafeMpi( MPI_Comm_compare( comm1, comm2, &result ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void
clique::mpi::ErrorHandlerSet
( mpi::Comm comm, mpi::ErrorHandler errorHandler )
{
#ifndef RELEASE
    PushCallStack("mpi::ErrorHandlerSet");
#endif
    SafeMpi( MPI_Errhandler_set( comm, errorHandler ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//---------------------------------//
// Cartesian communicator routines //
//---------------------------------//

void
clique::mpi::CartCreate
( mpi::Comm comm, int numDims, const int* dimensions, const int* periods, 
  int reorder, mpi::Comm& cartComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CartCreate");
#endif
    SafeMpi( 
        MPI_Cart_create
        ( comm, numDims, const_cast<int*>(dimensions), 
          const_cast<int*>(periods), reorder, &cartComm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::CartSub
( mpi::Comm comm, const int* remainingDims, mpi::Comm& subComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CartSub");
#endif
    SafeMpi( MPI_Cart_sub( comm, const_cast<int*>(remainingDims), &subComm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//--------------------//
// Group manipulation //
//--------------------//

int
clique::mpi::GroupRank( mpi::Group group )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupRank");
#endif
    int rank;
    SafeMpi( MPI_Group_rank( group, &rank ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

int
clique::mpi::GroupSize( mpi::Group group )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupSize");
#endif
    int size;
    SafeMpi( MPI_Group_size( group, &size ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return size;
}

void
clique::mpi::CommGroup( mpi::Comm comm, mpi::Group& group )
{
#ifndef RELEASE
    PushCallStack("mpi::CommGroup");
#endif
    SafeMpi( MPI_Comm_group( comm, &group ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::GroupIncl
( mpi::Group group, int n, const int* ranks, mpi::Group& subGroup )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupIncl");
#endif
    SafeMpi( MPI_Group_incl( group, n, const_cast<int*>(ranks), &subGroup ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::GroupDifference
( mpi::Group parent, mpi::Group subset, mpi::Group& complement )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupDifference");
#endif
    SafeMpi( MPI_Group_difference( parent, subset, &complement ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::GroupFree( mpi::Group& group )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupFree");
#endif
    SafeMpi( MPI_Group_free( &group ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::GroupTranslateRanks
( mpi::Group origGroup, int size, const int* origRanks, 
  mpi::Group newGroup,                  int* newRanks )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupTranslateRanks");
#endif
    SafeMpi( 
        MPI_Group_translate_ranks
        ( origGroup, size, const_cast<int*>(origRanks), newGroup, newRanks ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Wait until every process in comm reaches this statement
void
clique::mpi::Barrier( mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Barrier");
#endif
    SafeMpi( MPI_Barrier( comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Test for completion
int
clique::mpi::Test( mpi::Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::Test");
#endif
    MPI_Status status;
    int flag;
    SafeMpi( MPI_Test( &request, &flag, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return flag;
}

// Ensure that the request finishes before continuing
void
clique::mpi::Wait( mpi::Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::Wait");
#endif
    MPI_Status status;
    SafeMpi( MPI_Wait( &request, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Nonblocking test for message completion
void
clique::mpi::IProbe
( int source, int tag, MPI_Comm comm, int& flag, MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::IProbe");
#endif
    SafeMpi( MPI_Iprobe( source, tag, comm, &flag, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<>
int
clique::mpi::GetCount<char>( MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_CHAR, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int
clique::mpi::GetCount<int>( MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_INT, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int
clique::mpi::GetCount<float>( MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_FLOAT, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int
clique::mpi::GetCount<double>( MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_DOUBLE, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int
clique::mpi::GetCount<clique::scomplex>( MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_COMPLEX, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int
clique::mpi::GetCount<clique::dcomplex>( MPI_Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_DOUBLE_COMPLEX, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

void
clique::mpi::Send
( const char* buf, int count, int to, int tag, MPI_Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<char*>(buf), count, MPI_CHAR, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISend
( const char* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{ 
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<char*>(buf), count, MPI_CHAR, to, tag, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISSend
( const char* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<char*>(buf), count, MPI_CHAR, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Send
( const int* buf, int count, int to, int tag, MPI_Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<int*>(buf), count, MPI_INT, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISend
( const int* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{ 
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<int*>(buf), count, MPI_INT, to, tag, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISSend
( const int* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<int*>(buf), count, MPI_INT, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Send
( const float* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISend
( const float* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISSend
( const float* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Send
( const double* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISend
( const double* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISSend
( const double* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Send
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Send
        ( const_cast<scomplex*>(buf), 2*count, MPI_FLOAT, to, tag, comm )
    );
#else
    SafeMpi( 
        MPI_Send
        ( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISend
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Isend
        ( const_cast<scomplex*>(buf), 2*count, MPI_FLOAT, to, tag, comm,
          &request )
    );
#else
    SafeMpi( 
        MPI_Isend
        ( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISSend
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Issend
        ( const_cast<scomplex*>(buf), 2*count, MPI_FLOAT, to, tag, comm,
          &request )
    );
#else
    SafeMpi(
        MPI_Issend
        ( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Send
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Send
        ( const_cast<dcomplex*>(buf), 2*count, MPI_DOUBLE, to, tag, comm )
    );
#else
    SafeMpi( 
        MPI_Send
        ( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, tag, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISend
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Isend
        ( const_cast<dcomplex*>(buf), 2*count, MPI_DOUBLE, to, tag, comm,
          &request )
    );
#else
    SafeMpi( 
        MPI_Isend
        ( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ISSend
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Issend
        ( const_cast<dcomplex*>(buf), 2*count, MPI_DOUBLE, to, tag, comm,
          &request )
    );
#else
    SafeMpi(
        MPI_Issend
        ( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Recv
( char* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MPI_Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_CHAR, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::IRecv
( char* buf, int count, int from, int tag, MPI_Comm comm, 
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_CHAR, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Recv
( int* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MPI_Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_INT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::IRecv
( int* buf, int count, int from, int tag, MPI_Comm comm, 
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_INT, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Recv
( float* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MPI_Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_FLOAT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::IRecv
( float* buf, int count, int from, int tag, MPI_Comm comm, 
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_FLOAT, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Recv
( double* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MPI_Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_DOUBLE, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::IRecv
( double* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_DOUBLE, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Recv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MPI_Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Recv( buf, 2*count, MPI_FLOAT, from, tag, comm, &status ) );
#else
    SafeMpi( MPI_Recv( buf, count, MPI_COMPLEX, from, tag, comm, &status ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::IRecv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Irecv( buf, 2*count, MPI_FLOAT, from, tag, comm, &request ) );
#else
    SafeMpi( MPI_Irecv( buf, count, MPI_COMPLEX, from, tag, comm, &request ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Recv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MPI_Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Recv( buf, 2*count, MPI_DOUBLE, from, tag, comm, &status )
    );
#else
    SafeMpi( 
        MPI_Recv( buf, count, MPI_DOUBLE_COMPLEX, from, tag, comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::IRecv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Irecv( buf, 2*count, MPI_DOUBLE, from, tag, comm, &request )
    );
#else
    SafeMpi( 
        MPI_Irecv( buf, count, MPI_DOUBLE_COMPLEX, from, tag, comm, &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::SendRecv
( const char* sbuf, int sc, int to,   int stag,
        char* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    MPI_Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<char*>(sbuf), sc, MPI_CHAR, to, stag,
          rbuf, rc, MPI_CHAR, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::SendRecv
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    MPI_Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<int*>(sbuf), sc, MPI_INT, to, stag,
          rbuf, rc, MPI_INT, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::SendRecv
( const float* sbuf, int sc, int to,   int stag,
        float* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    MPI_Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT, to, stag,
          rbuf, rc, MPI_FLOAT, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::SendRecv
( const double* sbuf, int sc, int to,   int stag,
        double* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    MPI_Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE, to, stag,
          rbuf, rc, MPI_DOUBLE, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    MPI_Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Sendrecv
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT, to,   stag,
          rbuf,                        2*rc, MPI_FLOAT, from, rtag, 
          comm, &status )
    );
#else
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX, to,   stag,
          rbuf,                        rc, MPI_COMPLEX, from, rtag, 
          comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
        dcomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    MPI_Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Sendrecv
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE, to,   stag,
          rbuf,                        2*rc, MPI_DOUBLE, from, rtag, 
          comm, &status )
    );
#else
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX, to,   stag,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, from, rtag, 
          comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Broadcast
( char* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_CHAR, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Broadcast
( int* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Broadcast
( float* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Broadcast
( double* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Broadcast
( scomplex* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, MPI_FLOAT, root, comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, MPI_COMPLEX, root, comm ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Broadcast
( dcomplex* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, MPI_DOUBLE, root, comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, MPI_DOUBLE_COMPLEX, root, comm ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Gather
( const char* sbuf, int sc,
        char* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<char*>(sbuf), sc, MPI_CHAR,
          rbuf,                    rc, MPI_CHAR, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Gather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Gather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Gather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Gather
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, root, comm )
    );
#else
    SafeMpi( 
        MPI_Gather
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Gather
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, root, comm )
    );
#else
    SafeMpi( 
        MPI_Gather
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, root, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllGather
( const char* sbuf, int sc,
        char* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
    SafeMpi( 
        MPI_Allgather
        ( const_cast<char*>(sbuf), sc, MPI_CHAR, 
          rbuf,                    rc, MPI_CHAR, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_CHAR_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<int*>(sbuf), sizeof(int)*sc, MPI_CHAR, 
          rbuf,                   sizeof(int)*rc, MPI_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<int*>(sbuf), sc, MPI_INT, 
          rbuf,                   rc, MPI_INT, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_CHAR_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<float*>(sbuf), sizeof(float)*sc, MPI_CHAR, 
          rbuf,                     sizeof(float)*rc, MPI_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT, 
          rbuf,                     rc, MPI_FLOAT, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_CHAR_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<double*>(sbuf), sizeof(double)*sc, MPI_CHAR, 
          rbuf,                      sizeof(double)*rc, MPI_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_CHAR_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<scomplex*>(sbuf), 2*sizeof(float)*sc, MPI_CHAR, 
          rbuf,                        2*sizeof(float)*rc, MPI_CHAR, comm ) 
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Allgather
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, comm )
    );
 #else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_CHAR_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<dcomplex*>(sbuf), 2*sizeof(double)*sc, MPI_CHAR, 
          rbuf,                        2*sizeof(double)*rc, MPI_CHAR, comm ) 
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Allgather
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, comm )
    );
 #else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Scatter
( const char* sbuf, int sc,
        char* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<char*>(sbuf), sc, MPI_CHAR,
          rbuf,                    rc, MPI_CHAR, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Scatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Scatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Scatter
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, root, comm )
    );
#else
    SafeMpi( 
        MPI_Scatter
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Scatter
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, root, comm )
    );
#else
    SafeMpi( 
        MPI_Scatter
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const char* sbuf, int sc,
        char* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<char*>(sbuf), sc, MPI_CHAR,
          rbuf,                    rc, MPI_CHAR, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Alltoall
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Alltoall
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const char* sbuf, const int* scs, const int* sds, 
        char* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<char*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_CHAR,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_CHAR, 
          comm ) 
    ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const int* sbuf, const int* scs, const int* sds, 
        int* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<int*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_INT,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_INT, 
          comm ) 
    ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<float*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_FLOAT,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_FLOAT, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<double*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_DOUBLE,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    int p;
    MPI_Comm_size( comm, &p );
    std::vector<int> scsDoubled(p);
    std::vector<int> sdsDoubled(p);
    std::vector<int> rcsDoubled(p);
    std::vector<int> rdsDoubled(p);
    for( int i=0; i<p; ++i )
        scsDoubled[i] = 2*scs[i];
    for( int i=0; i<p; ++i )
        sdsDoubled[i] = 2*sds[i];
    for( int i=0; i<p; ++i )
        rcsDoubled[i] = 2*rcs[i];
    for( int i=0; i<p; ++i )
        rdsDoubled[i] = 2*rds[i];
    SafeMpi(
        MPI_Alltoallv
        ( const_cast<scomplex*>(sbuf),
                &scsDoubled[0], &sdsDoubled[0], MPI_FLOAT,
          rbuf, &rcsDoubled[0], &rdsDoubled[0], MPI_FLOAT, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<scomplex*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_COMPLEX,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_COMPLEX, 
          comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    int p;
    MPI_Comm_size( comm, &p );
    std::vector<int> scsDoubled(p);
    std::vector<int> sdsDoubled(p);
    std::vector<int> rcsDoubled(p);
    std::vector<int> rdsDoubled(p);
    for( int i=0; i<p; ++i )
        scsDoubled[i] = 2*scs[i];
    for( int i=0; i<p; ++i )
        sdsDoubled[i] = 2*sds[i];
    for( int i=0; i<p; ++i )
        rcsDoubled[i] = 2*rcs[i];
    for( int i=0; i<p; ++i )
        rdsDoubled[i] = 2*rds[i];
    SafeMpi(
        MPI_Alltoallv
        ( const_cast<dcomplex*>(sbuf),
                &scsDoubled[0], &sdsDoubled[0], MPI_DOUBLE,
          rbuf, &rcsDoubled[0], &rdsDoubled[0], MPI_DOUBLE, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<dcomplex*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_DOUBLE_COMPLEX,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE_COMPLEX, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Reduce
( const char* sbuf, char* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<char*>(sbuf), rbuf, count, MPI_CHAR, op, root, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Reduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<int*>(sbuf), rbuf, count, MPI_INT, op, root, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Reduce
( const float* sbuf, float* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<float*>(sbuf), rbuf, count, MPI_FLOAT, op, root, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Reduce
( const double* sbuf, double* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<double*>(sbuf), rbuf, count, MPI_DOUBLE, op, root, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == MPI_SUM )
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, 2*count, MPI_FLOAT, op, root, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, count, MPI_COMPLEX, op, root, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Reduce
            ( const_cast<scomplex*>(sbuf), 
              rbuf, count, MPI_COMPLEX, op, root, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == MPI_SUM )
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, 2*count, MPI_DOUBLE, op, root, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, count, MPI_DOUBLE_COMPLEX, op, root, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Reduce
            ( const_cast<dcomplex*>(sbuf), 
              rbuf, count, MPI_DOUBLE_COMPLEX, op, root, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllReduce
( const char* sbuf, char* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<char*>(sbuf), rbuf, count, MPI_CHAR, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllReduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<int*>(sbuf), rbuf, count, MPI_INT, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllReduce
( const float* sbuf, float* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<float*>(sbuf), rbuf, count, MPI_FLOAT, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllReduce
( const double* sbuf, double* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<double*>(sbuf), rbuf, count, MPI_DOUBLE, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == MPI_SUM )
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, 2*count, MPI_FLOAT, op, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, count, MPI_COMPLEX, op, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<scomplex*>(sbuf), 
              rbuf, count, MPI_COMPLEX, op, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == MPI_SUM )
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, 2*count, MPI_DOUBLE, op, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, count, MPI_DOUBLE_COMPLEX, op, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<dcomplex*>(sbuf), 
              rbuf, count, MPI_DOUBLE_COMPLEX, op, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ReduceScatter
( const char* sbuf, char* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<char*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_CHAR, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<int*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_INT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<float*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_FLOAT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<double*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_DOUBLE, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, MPI_Op op, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    if( op == MPI_SUM )
    {
        int p;
        MPI_Comm_size( comm, &p );
        std::vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<scomplex*>(sbuf),
              rbuf, &rcsDoubled[0], MPI_FLOAT, op, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<scomplex*>(sbuf),
              rbuf, const_cast<int*>(rcs), MPI_COMPLEX, op, comm )
        );
    }
#else
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<scomplex*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_COMPLEX, op, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::mpi::ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, MPI_Op op, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    if( op == MPI_SUM )
    {
        int p;
        MPI_Comm_size( comm, &p );
        std::vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<dcomplex*>(sbuf),
              rbuf, &rcsDoubled[0], MPI_DOUBLE, op, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<dcomplex*>(sbuf),
              rbuf, const_cast<int*>(rcs), MPI_DOUBLE_COMPLEX, op, comm )
        );
    }
#else
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<dcomplex*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_DOUBLE_COMPLEX, op, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

