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

template<typename F>
void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize, 
  const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerMultiplyNormal");
    if( L.Height() != L.Width() || L.Height() != X.Height() || 
        L.Height() < supernodeSize )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  supernodeSize ~ " << supernodeSize << "\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    // Matrix views
    Matrix<F>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

    // Start the algorithm
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XT.Height() < supernodeSize )
    {
        const int blocksize = std::min(Blocksize(),supernodeSize-XT.Height());
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22, blocksize );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2, blocksize );

        //--------------------------------------------------------------------//
        // HERE
        throw std::logic_error("This routine is not yet finished");
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template<typename F>
void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset,
  int supernodeSize, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    clique::PushCallStack("numeric::LocalFrontLowerMultiplyTranspose");
    if( L.Height() != L.Width() || L.Height() != X.Height() || 
        L.Height() < supernodeSize )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  supernodeSize ~ " << supernodeSize << "\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transposed");
#endif
    // Matrix views
    Matrix<F>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

    Matrix<F> XT,  X0,
              XB,  X1,
                   X2;

    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, L.Height()-supernodeSize );
    PartitionUp
    ( X, XT,
         XB, X.Height()-supernodeSize );

    throw std::logic_error("This routine is not yet finished");

    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,   L00, L01, /**/ L02,
               /**/        L10, L11, /**/ L12,
         /*************/  /******************/
          LBL, /**/ LBR,   L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        //--------------------------------------------------------------------//
        // HERE
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    clique::PopCallStack();
#endif
}

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<float>& L, Matrix<float>& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<float>& L, Matrix<float>& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<double>& L, Matrix<double>& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<double>& L, Matrix<double>& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<float> >& L, Matrix<std::complex<float> >& X );

template void clique::numeric::LocalFrontLowerMultiplyNormal
( Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
template void clique::numeric::LocalFrontLowerMultiplyTranspose
( Orientation orientation, Diagonal diag, int diagOffset, int supernodeSize,
  const Matrix<std::complex<double> >& L, Matrix<std::complex<double> >& X );
