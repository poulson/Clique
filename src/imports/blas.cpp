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

//----------------------------------------------------------------------------//
// Level 1 BLAS                                                               //
//----------------------------------------------------------------------------//
void
clique::blas::Axpy
( int n, int alpha, const int* x, int incx, int* y, int incy )
{
    for( int i=0; i<n; ++i )
        y[i*incy] += alpha*x[i*incx];
}

void
clique::blas::Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy )
{ BLAS(saxpy)( &n, &alpha, x, &incx, y, &incy ); }

void
clique::blas::Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy )
{ BLAS(daxpy)( &n, &alpha, x, &incx, y, &incy ); }

void
clique::blas::Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy )
{ BLAS(caxpy)( &n, &alpha, x, &incx, y, &incy ); }

void
clique::blas::Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy )
{ BLAS(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }

float
clique::blas::Dot
( int n, const float* x, int incx, const float* y, int incy )
{ return BLAS(sdot)( &n, x, &incx, y, &incy ); }

double
clique::blas::Dot
( int n, const double* x, int incx, const double* y, int incy )
{ return BLAS(ddot)( &n, x, &incx, y, &incy ); }

clique::scomplex
clique::blas::Dot
( int n, const clique::scomplex* x, int incx,
         const clique::scomplex* y, int incy )
{ 
    clique::scomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += std::conj(x[i*incx])*y[i*incy];
    return alpha;
}

clique::dcomplex
clique::blas::Dot
( int n, const clique::dcomplex* x, int incx,
         const clique::dcomplex* y, int incy )
{
    clique::dcomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += std::conj(x[i*incx])*y[i*incy];
    return alpha;
}

float
clique::blas::Dotc
( int n, const float* x, int incx, const float* y, int incy )
{ return BLAS(sdot)( &n, x, &incx, y, &incy ); }

double
clique::blas::Dotc
( int n, const double* x, int incx, const double* y, int incy )
{ return BLAS(ddot)( &n, x, &incx, y, &incy ); }

clique::scomplex
clique::blas::Dotc
( int n, const clique::scomplex* x, int incx,
         const clique::scomplex* y, int incy )
{ 
    clique::scomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += std::conj(x[i*incx])*y[i*incy];
    return alpha;
}

clique::dcomplex
clique::blas::Dotc
( int n, const clique::dcomplex* x, int incx,
         const clique::dcomplex* y, int incy )
{ 
    clique::dcomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += std::conj(x[i*incx])*y[i*incy];
    return alpha;
}

float
clique::blas::Dotu
( int n, const float* x, int incx, const float* y, int incy )
{ return BLAS(sdot)( &n, x, &incx, y, &incy ); }

double
clique::blas::Dotu
( int n, const double* x, int incx, const double* y, int incy )
{ return BLAS(ddot)( &n, x, &incx, y, &incy ); }

clique::scomplex
clique::blas::Dotu
( int n, const clique::scomplex* x, int incx,
         const clique::scomplex* y, int incy )
{
    clique::scomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}

clique::dcomplex
clique::blas::Dotu
( int n, const clique::dcomplex* x, int incx,
         const clique::dcomplex* y, int incy )
{
    clique::dcomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}

float
clique::blas::Nrm2
( int n, const float* x, int incx )
{ return BLAS(snrm2)( &n, x, &incx ); }

double
clique::blas::Nrm2
( int n, const double* x, int incx )
{ return BLAS(dnrm2)( &n, x, &incx ); }

float
clique::blas::Nrm2
( int n, const scomplex* x, int incx )
{ return BLAS(scnrm2)( &n, x, &incx ); }

double
clique::blas::Nrm2
( int n, const dcomplex* x, int incx )
{ return BLAS(dznrm2)( &n, x, &incx ); }

void
clique::blas::Scal
( int n, float alpha, float* x, int incx )
{ BLAS(sscal)( &n, &alpha, x, &incx ); }

void
clique::blas::Scal
( int n, double alpha, double* x, int incx )
{ BLAS(dscal)( &n, &alpha, x, &incx ); }

void
clique::blas::Scal
( int n, scomplex alpha, scomplex* x, int incx )
{ BLAS(cscal)( &n, &alpha, x, &incx ); }

void
clique::blas::Scal
( int n, dcomplex alpha, dcomplex* x, int incx )
{ BLAS(zscal)( &n, &alpha, x, &incx ); }

//----------------------------------------------------------------------------//
// Level 2 BLAS                                                               //
//----------------------------------------------------------------------------//
void
clique::blas::Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void
clique::blas::Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void
clique::blas::Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ BLAS(cgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ BLAS(zgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda  )
{ BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ BLAS(chemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ BLAS(zhemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda )
{ BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void
clique::blas::Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void
clique::blas::Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{ BLAS(cher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void
clique::blas::Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{ BLAS(zher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void
clique::blas::Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Symv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void
clique::blas::Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    LAPACK(csymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void
clique::blas::Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    LAPACK(zsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void
clique::blas::Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda  )
{ BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void
clique::blas::Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void
clique::blas::Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    LAPACK(csyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void
clique::blas::Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    LAPACK(zsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void
clique::blas::Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void
clique::blas::Syr2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{
    // csyr2 doesn't exist, so we route through csyr2k. However, csyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const int k = 1;
    const scomplex beta = 1.;
    BLAS(csyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

void
clique::blas::Syr2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{
    // zsyr2 doesn't exist, so we route through zsyr2k. However, zsyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const int k = 1;
    const dcomplex beta = 1.;
    BLAS(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

void
clique::blas::Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ BLAS(strmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ BLAS(dtrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ BLAS(ctrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ BLAS(ztrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ BLAS(strsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ BLAS(dtrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ BLAS(ctrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void
clique::blas::Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ BLAS(ztrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

//----------------------------------------------------------------------------//
// Level 3 BLAS                                                               //
//----------------------------------------------------------------------------//
void
clique::blas::Gemm
( char transA, char transB, int m, int n, int k, 
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    BLAS(sgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Gemm
( char transA, char transB,
  int m, int n, int k, 
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    BLAS(dgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Gemm
( char transA, char transB, int m, int n, int k, 
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(cgemm)( &transA, &transB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Gemm
( char transA, char transB, int m, int n, int k, 
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zgemm)( &transA, &transB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    BLAS(ssymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    BLAS(dsymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(chemm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zhemm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(ssyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(dsyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(cher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(ssyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void
clique::blas::Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(dsyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void
clique::blas::Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ BLAS(cherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void
clique::blas::Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ BLAS(zherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void
clique::blas::Hetrmm( char uplo, int n, float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("blas::Hetrmm");
#endif
    int info;
    LAPACK(slauum)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream os;
        os << "slauum returned with info=" << info;
        throw std::logic_error( os.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::blas::Hetrmm( char uplo, int n, double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("blas::Hetrmm");
#endif
    int info;
    LAPACK(dlauum)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream os;
        os << "dlauum returned with info=" << info;
        throw std::logic_error( os.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::blas::Hetrmm( char uplo, int n, scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("blas::Hetrmm");
#endif
    int info;
    LAPACK(clauum)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream os;
        os << "clauum returned with info=" << info;
        throw std::logic_error( os.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::blas::Hetrmm( char uplo, int n, dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("blas::Hetrmm");
#endif
    int info;
    LAPACK(zlauum)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream os;
        os << "zlauum returned with info=" << info;
        throw std::logic_error( os.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
clique::blas::Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    BLAS(ssymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    BLAS(dsymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(csymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zsymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    BLAS(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    BLAS(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void
clique::blas::Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{ BLAS(ssyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void
clique::blas::Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{ BLAS(dsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void
clique::blas::Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ BLAS(csyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void
clique::blas::Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ BLAS(zsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void
clique::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    BLAS(strmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void
clique::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    BLAS(dtrmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void
clique::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    BLAS(ctrmm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void
clique::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    BLAS(ztrmm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void
clique::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(strsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

void
clique::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(dtrsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

void
clique::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    BLAS(ctrsm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

void
clique::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    BLAS(ztrsm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

