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

#ifndef WITHOUT_PMRRR

namespace elem {

//
// Invert the sufficiently large eigenvalues of A.
//

template<typename F>
inline void
HermitianPseudoinverse
( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianPseudoinverse");
#endif
    typedef typename Base<F>::type R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Compute the two-norm of A as the maximum absolute value of its
    // eigenvalues
    R maxLocalAbsEig = 0;
    const int numLocalEigs = w.LocalHeight();
    for( int iLocal=0; iLocal<numLocalEigs; ++iLocal )
    {
        const R omega = w.GetLocal(iLocal,0);
        maxLocalAbsEig = std::max(maxLocalAbsEig,Abs(omega));
    }
    R twoNorm;
    mpi::AllReduce( &maxLocalAbsEig, &twoNorm, 1, mpi::MAX, g.VCComm() );

    // Set the tolerance equal to n ||A||_2 eps, and invert values above it
    const int n = A.Height();
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = n*twoNorm*eps;
    for( int iLocal=0; iLocal<numLocalEigs; ++iLocal )
    {
        const R omega = w.GetLocal(iLocal,0);
        if( Abs(omega) < tolerance )
            w.SetLocal(iLocal,0,0);
        else
            w.SetLocal(iLocal,0,1/omega);
    }

    // Form the pseudoinverse
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR