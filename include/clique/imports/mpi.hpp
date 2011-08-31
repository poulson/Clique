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
#ifndef CLIQUE_IMPORTS_MPI_HPP
#define CLIQUE_IMPORTS_MPI_HPP 1

namespace clique {
namespace mpi {

// Datatype definitions
typedef MPI_Comm Comm;
typedef MPI_Datatype Datatype;
typedef MPI_Errhandler ErrorHandler;
typedef MPI_Group Group;
typedef MPI_Op Op;
typedef MPI_Request Request;
typedef MPI_Status Status;
typedef MPI_User_function UserFunction;

// Standard constants
const int ANY_SOURCE = MPI_ANY_SOURCE;
const int ANY_TAG = MPI_ANY_TAG;
const int THREAD_SINGLE = MPI_THREAD_SINGLE;
const int THREAD_FUNNELED = MPI_THREAD_FUNNELED;
const int THREAD_SERIALIZED = MPI_THREAD_SERIALIZED;
const int THREAD_MULTIPLE = MPI_THREAD_MULTIPLE;
const int UNDEFINED = MPI_UNDEFINED;
const Comm COMM_WORLD = MPI_COMM_WORLD;
const ErrorHandler ERRORS_RETURN = MPI_ERRORS_RETURN;
const ErrorHandler ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL;
const Group GROUP_EMPTY = MPI_GROUP_EMPTY;
const Op MAX = MPI_MAX;
const Op SUM = MPI_SUM;

// Added constant(s)
const int MIN_COLL_MSG = 1; // minimum message size for collectives

//----------------------------------------------------------------------------//
// Routines                                                                   //
//----------------------------------------------------------------------------//

// MPI environmental routines
void Init( int& argc, char**& argv );
int InitThread( int& argc, char**& argv, int required );
void Finalize();
int Initialized();
int Finalized();
double Time();
void OpCreate( UserFunction* func, int commutes, Op& op );
void OpFree( Op& op );

// Communicator manipulation
int CommRank( Comm comm );
int CommSize( Comm comm );
void CommCreate( Comm parentComm, Group subsetGroup, Comm& subsetComm );
void CommDup( Comm original, Comm& duplicate );
void CommSplit( Comm comm, int color, int key, Comm& newComm );
void CommFree( Comm& comm );
bool CongruentComms( Comm comm1, Comm comm2 );
void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler );

// Cartesian communicator routines
void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  int reorder, Comm& cartComm );
void CartSub
( Comm comm, const int* remainingDims, Comm& subComm );

// Group manipulation
int GroupRank( Group group );
int GroupSize( Group group );
void CommGroup( Comm comm, Group& group );
void GroupIncl( Group group, int n, const int* ranks, Group& subGroup );
void GroupDifference( Group parent, Group subset, Group& complement );
void GroupFree( Group& group );
void GroupTranslateRanks
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks );

// Wait until every process in comm reaches this statement
void Barrier( Comm comm );

// Ensure that the request finishes before continuing
void Wait( Request& request );

// Test for completion
int Test( Request& request );

// Nonblocking test for message completion
void IProbe
( int source, int tag, Comm comm, int& flag, Status& status );

template<typename T>
int GetCount( Status& status );

template<> int GetCount<byte>( Status& status );
template<> int GetCount<int>( Status& status );
template<> int GetCount<float>( Status& status );
template<> int GetCount<double>( Status& status );
template<> int GetCount<scomplex>( Status& status );
template<> int GetCount<dcomplex>( Status& status );

void
Send
( const byte* buf, int count, int to, int tag, Comm comm );
void
ISend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
void
ISSend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request );

void
Send
( const int* buf, int count, int to, int tag, Comm comm );
void
ISend
( const int* buf, int count, int to, int tag, Comm comm, Request& request );
void
ISSend
( const int* buf, int count, int to, int tag, Comm comm, Request& request );

void
Send
( const float* buf, int count, int to, int tag, Comm comm );
void
ISend
( const float* buf, int count, int to, int tag, Comm comm, Request& request );
void
ISSend
( const float* buf, int count, int to, int tag, Comm comm, Request& request );

void 
Send
( const double* buf, int count, int to, int tag, Comm comm );
void
ISend
( const double* buf, int count, int to, int tag, Comm comm, Request& request );
void
ISSend
( const double* buf, int count, int to, int tag, Comm comm, Request& request );

void
Send
( const scomplex* buf, int count, int to, int tag, Comm comm );
void
ISend
( const scomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
void
ISSend
( const scomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );

void
Send
( const dcomplex* buf, int count, int to, int tag, Comm comm );
void
ISend
( const dcomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
void
ISSend
( const dcomplex* buf, int count, int to, int tag, Comm comm, 
 Request& request );
 
void
Recv
( byte* buf, int count, int from, int tag, Comm comm );
void
IRecv
( byte* buf, int count, int from, int tag, Comm comm, Request& request );
       
void
Recv
( int* buf, int count, int from, int tag, Comm comm );
void
IRecv
( int* buf, int count, int from, int tag, Comm comm, Request& request );

void
Recv
( float* buf, int count, int from, int tag, Comm comm );
void
IRecv
( float* buf, int count, int from, int tag, Comm comm, Request& request );

void
Recv
( double* buf, int count, int from, int tag, Comm comm );
void
IRecv
( double* buf, int count, int from, int tag, Comm comm, Request& request );

void
Recv
( scomplex* buf, int count, int from, int tag, Comm comm );
void
IRecv
( scomplex* buf, int count, int from, int tag, Comm comm, Request& request );

void
Recv
( dcomplex* buf, int count, int from, int tag, Comm comm );
void
IRecv
( dcomplex* buf, int count, int from, int tag, Comm comm, Request& request );

void
SendRecv
( const byte* sbuf, int sc, int to,   int stag,
        byte* rbuf, int rc, int from, int rtag, Comm comm );

void
SendRecv
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, Comm comm );

void
SendRecv
( const float* sbuf, int sc, int to,   int stag,
        float* rbuf, int rc, int from, int rtag, Comm comm );

void
SendRecv
( const double* sbuf, int sc, int to,   int stag,
        double* rbuf, int rc, int from, int rtag, Comm comm );

void
SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, Comm comm );

void
SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
        dcomplex* rbuf, int rc, int from, int rtag, Comm comm );

void
Broadcast( byte* buf, int count, int root, Comm comm );

void
Broadcast( int* buf, int count, int root, Comm comm );

void
Broadcast( float* buf, int count, int root, Comm comm );

void
Broadcast( double* buf, int count, int root, Comm comm );

void
Broadcast
( scomplex* buf, int count, int root, Comm comm );

void
Broadcast
( dcomplex* buf, int count, int root, Comm comm );

void
Gather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm );

void
Gather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm );

void
Gather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm );

void
Gather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm );

void 
Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm );

void
Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm );
 
void
AllGather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm );
   
void
AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm );

void
AllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm );

void
AllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm );

void
AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm );

void
AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm );

void
Scatter
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm );

void
Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm );

void
Scatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm );

void
Scatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm );

void
Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm );

void
Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm );
 
void
AllToAll
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm );
   
void
AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm );

void
AllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm );

void
AllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm );

void
AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm );

void
AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm );

void
AllToAll
( const byte* sbuf, const int* scs, const int* sds,
        byte* rbuf, const int* rcs, const int* rds, Comm comm );

void
AllToAll
( const int* sbuf, const int* scs, const int* sds,
        int* rbuf, const int* rcs, const int* rds, Comm comm );

void
AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, Comm comm );

void
AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, Comm comm );

void
AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, Comm comm );

void
AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, Comm comm );

void
Reduce
( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );

void
Reduce
( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );

void
Reduce
( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );

void
Reduce
( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );

void
Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, int root, Comm comm );

void
Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, int root, Comm comm );
    
void
AllReduce
( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );

void
AllReduce
( const int* sbuf, int* rbuf, int count, Op op, Comm comm );

void
AllReduce
( const float* sbuf, float* rbuf, int count, Op op, Comm comm );

void
AllReduce
( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
void
AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, Comm comm );

void
AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, Comm comm );
        
void
ReduceScatter
( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm );

void
ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm );

void
ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm );

void
ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm );
void
ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, Op op, Comm comm );

void
ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, Op op, Comm comm );

} // mpi
} // clique

#endif /* CLIQUE_IMPORTS_MPI_HPP */

