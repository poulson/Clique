/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

const int n = 300;

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and 
    // does so if necessary. The full routine is elem::Initialize.
    Initialize( argc, argv );

    // Extract our MPI rank
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    // Surround the Elemental calls with try/catch statements in order to 
    // safely handle any exceptions that were thrown during execution.
    try 
    {
        // Create a 2d process grid from a communicator. In our case, it is
        // MPI_COMM_WORLD. There is another constructor that allows you to 
        // specify the grid dimensions, Grid g( comm, r, c ), which creates an 
        // r x c grid.
        Grid g( comm );
    
        // Create an n x n complex distributed matrix, 
        // We distribute the matrix using grid 'g'.
        //
        // There are quite a few available constructors, including ones that 
        // allow you to pass in your own local buffer and to specify the 
        // distribution alignments (i.e., which process row and column owns the
        // top-left element)
        DistMatrix<C> H( n, n, g );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
        // the global matrix is Hermitian. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        //
        const int colShift = H.ColShift(); // first row we own
        const int rowShift = H.RowShift(); // first col we own
        const int colStride = H.ColStride();
        const int rowStride = H.RowStride();
        const int localHeight = H.LocalHeight();
        const int localWidth = H.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                // Our process owns the rows colShift:colStride:n,
                //           and the columns rowShift:rowStride:n
                const int i = colShift + iLocal*colStride;
                const int j = rowShift + jLocal*rowStride;
                H.SetLocal( iLocal, jLocal, C(i+j,i-j) );
            }
        }

        // Make a backup of H before we overwrite it within the eigensolver
        DistMatrix<C> HCopy( H );

        // Call the eigensolver. We first create an empty complex eigenvector 
        // matrix, X[MC,MR], and an eigenvalue column vector, w[VR,* ]
        //
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        DistMatrix<R,VR,STAR> w( g );
        DistMatrix<C> X( g );
        HermitianEig( LOWER, H, w, X ); // only use lower half of H

        // Optional: sort the eigenpairs
        SortEig( w, X );

        // Check the residual, || H X - Omega X ||_F
        const R frobH = HermitianNorm( LOWER, HCopy, FROBENIUS_NORM );
        DistMatrix<C> E( X );
        DiagonalScale( RIGHT, NORMAL, w, E );
        Hemm( LEFT, LOWER, (C)-1, HCopy, X, (C)1, E );
        const R frobResid = HermitianNorm( LOWER, E, FROBENIUS_NORM );

        // Check the orthogonality of X
        Identity( n, n, E );
        Herk( LOWER, NORMAL, (C)-1, X, (C)1, E );
        const R frobOrthog = HermitianNorm( LOWER, E, FROBENIUS_NORM );

        if( g.Rank() == 0 )
        {
            std::cout << "|| H ||_F = " << frobH << "\n"
                      << "|| H X - X Omega ||_F / || A ||_F = " 
                      << frobResid / frobH << "\n"
                      << "|| X X^H - I ||_F = " << frobOrthog / frobH
                      << "\n" << std::endl;
        }
    }
    catch( exception& e )
    {
        cerr << "Process " << commRank << " caught exception with message: "
             << e.what() << endl;
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
