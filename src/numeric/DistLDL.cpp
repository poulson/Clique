/*
   Clique: a scalable implementation of the multifrontal algorithm

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
using namespace elemental;

template<typename F> // F represents a real or complex field
void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S, // can't be const due to map...
  const numeric::LocalSymmFact<F>& localL,
        numeric::DistSymmFact<F>& distL )
{
#ifndef RELEASE
    PushCallStack("numeric::DistLDL");
    if( orientation == NORMAL )
        throw std::logic_error("LDL must be (conjugate-)transposed");
#endif
    const int numSupernodes = S.supernodes.size();
    if( numSupernodes == 0 )
        return;

    // The bottom front is already computed, so just view it
    distL.fronts[0].LocalMatrix().LockedView( localL.fronts.back() );

    // Perform the distributed portion of the factorization
    for( unsigned k=1; k<numSupernodes; ++k )
    {
        const symbolic::DistSymmFactSupernode& symbSN = S.supernodes[k];
        DistMatrix<F,MC,MR>& front = distL.fronts[k];

        const Grid& g = front.Grid();
        mpi::Comm comm = g.VCComm();
        const unsigned commRank = mpi::CommRank( comm );
        const unsigned commSize = mpi::CommSize( comm );

#ifndef RELEASE
        if( front.Height() != symbSN.size+symbSN.lowerStruct.size() ||
            front.Width()  != symbSN.size+symbSN.lowerStruct.size() )
            throw std::logic_error("Front was not the proper size");
#endif

        // Pack our child's updates
        const DistMatrix<F,MC,MR>& childFront = distL.fronts[k-1];
        const bool isLeftChild = ( commRank < commSize/2 );
        // TODO: Use relative indices to compute the maximum size that each
        //       process will contribute, then pack our child's updates for
        //       the AllToAll.

        // AllToAll to send and receive the child updates
        // TODO

        // Unpack the child udpates (with an Axpy)
        // TODO

        DistSupernodeLDL( orientation, front, symbSN.size );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<float>& localL,
        numeric::DistSymmFact<float>& distL );

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<double>& localL,
        numeric::DistSymmFact<double>& distL );

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<float> >& localL,
        numeric::DistSymmFact<std::complex<float> >& distL );

template void clique::numeric::DistLDL
( Orientation orientation,
        symbolic::DistSymmFact& S,
  const numeric::LocalSymmFact<std::complex<double> >& localL,
        numeric::DistSymmFact<std::complex<double> >& distL );

