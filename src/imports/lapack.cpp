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

// Machine constants
template<> float 
clique::lapack::MachineEpsilon<float>()
{
    const char cmach = 'E';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineEpsilon<double>()
{
    const char cmach = 'E';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachineSafeMin<float>()
{
    const char cmach = 'S';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineSafeMin<double>()
{
    const char cmach = 'S';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachineBase<float>()
{
    const char cmach = 'B';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineBase<double>()
{
    const char cmach = 'B';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachinePrecision<float>()
{
    const char cmach = 'P';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachinePrecision<double>()
{
    const char cmach = 'P';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachineUnderflowExponent<float>()
{
    const char cmach = 'M';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineUnderflowExponent<double>()
{
    const char cmach = 'M';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachineUnderflowThreshold<float>()
{
    const char cmach = 'U';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineUnderflowThreshold<double>()
{
    const char cmach = 'U';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachineOverflowExponent<float>()
{
    const char cmach = 'L';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineOverflowExponent<double>()
{
    const char cmach = 'L';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
clique::lapack::MachineOverflowThreshold<float>()
{
    const char cmach = 'O';
    return LAPACK(slamch)( &cmach );
}

template<> double 
clique::lapack::MachineOverflowThreshold<double>()
{
    const char cmach = 'O';
    return LAPACK(dlamch)( &cmach );
}

void
clique::lapack::Chol
( char uplo, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    int info;
    LAPACK(spotrf)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "spotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::Chol
( char uplo, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    int info;
    LAPACK(dpotrf)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::Chol
( char uplo, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    int info;
    LAPACK(cpotrf)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::Chol
( char uplo, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Chol");
#endif
    int info;
    LAPACK(zpotrf)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::LU
( int m, int n, float* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(sgetrf)( &m, &n, A, &lda, p, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::LU
( int m, int n, double* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(dgetrf)( &m, &n, A, &lda, p, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::LU
( int m, int n, scomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(cgetrf)( &m, &n, A, &lda, p, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::lapack::LU
( int m, int n, dcomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(zgetrf)( &m, &n, A, &lda, p, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

float
clique::lapack::SafeNorm
( float alpha, float beta )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    float gamma = LAPACK(slapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

double
clique::lapack::SafeNorm
( double alpha, double beta )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    double gamma = LAPACK(dlapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

float
clique::lapack::SafeNorm
( float alpha, float beta, float gamma )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    float delta = LAPACK(slapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

double
clique::lapack::SafeNorm
( double alpha, double beta, double gamma )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    double delta = LAPACK(dlapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

