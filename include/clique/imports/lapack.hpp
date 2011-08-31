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
#ifndef CLIQUE_LAPACK_HPP
#define CLIQUE_LAPACK_HPP 1

namespace clique {
namespace lapack {

// Relative machine precision
template<typename R> R MachineEpsilon();
template<> float MachineEpsilon<float>();
template<> double MachineEpsilon<double>();

// Minimum number which can be inverted without overflow
template<typename R> R MachineSafeMin();
template<> float MachineSafeMin<float>();
template<> double MachineSafeMin<double>();

// Base of the machine, where the number is represented as 
//   (mantissa) x (base)^(exponent)
template<typename R> R MachineBase();
template<> float MachineBase<float>();
template<> double MachineBase<double>();

// Return the relative machine precision multiplied by the base
template<typename R> R MachinePrecision();
template<> float MachinePrecision<float>();
template<> double MachinePrecision<double>();

// Return the minimum exponent before (gradual) underflow occurs
template<typename R> R MachineUnderflowExponent();
template<> float MachineUnderflowExponent<float>();
template<> double MachineUnderflowExponent<double>();

// Return the underflow threshold: (base)^((underflow exponent)-1)
template<typename R> R MachineUnderflowThreshold();
template<> float MachineUnderflowThreshold<float>();
template<> double MachineUnderflowThreshold<double>();

// Return the largest exponent before overflow
template<typename R> R MachineOverflowExponent();
template<> float MachineOverflowExponent<float>();
template<> double MachineOverflowExponent<double>();

// Return the overflow threshold: (1-(rel. prec.)) * (base)^(overflow exponent)
template<typename R> R MachineOverflowThreshold();
template<> float MachineOverflowThreshold<float>();
template<> double MachineOverflowThreshold<double>();

void Chol( char uplo, int n, const float* A, int lda );
void Chol( char uplo, int n, const double* A, int lda );
void Chol( char uplo, int n, const scomplex* A, int lda );
void Chol( char uplo, int n, const dcomplex* A, int lda );

void LU( int m, int n, float* A, int lda, int* p );
void LU( int m, int n, double* A, int lda, int* p );
void LU( int m, int n, scomplex* A, int lda, int* p );
void LU( int m, int n, dcomplex* A, int lda, int* p );

float SafeNorm( float alpha, float beta );
double SafeNorm( double alpha, double beta );
float SafeNorm( float alpha, float beta, float gamma );
double SafeNorm( double alpha, double beta, double gamma );

} // lapack
} // clique

extern "C" {

// Machine constants
float LAPACK(slamch)( const char* cmach );
double LAPACK(dlamch)( const char* cmach );

// LU factorization
void LAPACK(sgetrf)
( const int* m, const int* n, 
  float* A, const int* lda, int* p, int* info );

void LAPACK(dgetrf)
( const int* m, const int* n, 
  double* A, const int* lda, int* p, int* info );

void LAPACK(cgetrf)
( const int* m, const int* n, 
  clique::scomplex* A, const int* lda, int* p, int* info );

void LAPACK(zgetrf)
( const int* m, const int* n, 
  clique::dcomplex* A, const int* lda, int* p, int* info );

// Safe norms
float LAPACK(slapy2)
( const float* alpha, const float* beta );

double LAPACK(dlapy2)
( const double* alpha, const double* beta );

float LAPACK(slapy3)
( const float* alpha, const float* beta, const float* gamma );

double LAPACK(dlapy3)
( const double* alpha, const double* beta, const double* gamma );

// Cholesky factorization
void LAPACK(spotrf)
( const char* uplo, const int* n, const float* A, const int* lda,
  int* info );

void LAPACK(dpotrf)
( const char* uplo, const int* n, const double* A, const int* lda,
  int* info );
    
void LAPACK(cpotrf)
( const char* uplo, const int* n, const clique::scomplex* A, 
  const int* lda, int* info );
    
void LAPACK(zpotrf)
( const char* uplo, const int* n, const clique::dcomplex* A, 
  const int* lda, int* info );

} // extern "C"

#endif /* CLIQUE_LAPACK_HPP */
