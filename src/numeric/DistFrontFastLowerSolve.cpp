/*
   Modification of include/elemental/basic/level3/Trsm/TrsmLLN.hpp 
   from Elemental.
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
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

namespace cliq {

template<typename F>
void numeric::DistFrontFastLowerForwardSolve
( Diagonal diag, DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontFastLowerForwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
#endif
    const Grid& g = L.Grid();
    const int commRank = g.VCRank();
    const int commSize = g.Size();
    if( commSize == 1 )
    {
        numeric::LocalFrontLowerForwardSolve
        ( diag, L.LockedLocalMatrix(), X.LocalMatrix() );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    PartitionDown
    ( L, LT, 
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    const int localTopHeight = LT.LocalHeight();
    std::vector<F> localDiag;
    if( diag == UNIT )
    {
        // Extract the diagonal of the top triangle and replace it with ones
        localDiag.resize( localTopHeight );
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            localDiag[iLocal] = LTBuffer[iLocal+i*LTLDim];
            LTBuffer[iLocal+i*LTLDim] = 1;
        }
    }

    // Get a copy of all of XT
    DistMatrix<F,STAR,STAR> XT_STAR_STAR( XT );

    // XT := LT XT
    elem::internal::LocalGemm
    ( NORMAL, NORMAL, (F)1, LT, XT_STAR_STAR, (F)0, XT );

    if( diag == UNIT )
    {
        // Put the diagonal back
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            LTBuffer[iLocal+i*LTLDim] = localDiag[iLocal];
        }
    }

    if( LB.Height() != 0 )
    {
        // Gather all of XT again
        XT_STAR_STAR = XT;

        // XB := XB - LB XT
        elem::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, LB, XT_STAR_STAR, (F)1, XB );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void numeric::DistFrontFastLowerForwardSolve
( Diagonal diag, DistMatrix<F,MC,MR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontFastLowerForwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = L.Grid();
    const int commRank = g.VCRank();
    const int commSize = g.Size();
    if( commSize == 1 )
    {
        numeric::LocalFrontLowerForwardSolve
        ( diag, L.LockedLocalMatrix(), X.LocalMatrix() );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F,MC,MR> LT(g),
                        LB(g);
    PartitionDown
    ( L, LT, 
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    DistMatrix<F,MD,STAR> dTOrig(g), dTReplacement(g);
    if( diag == UNIT )
    {
        // Extract the diagonal of the top triangle and replace it with ones
        LT.GetDiagonal( dTOrig );
        const int localHeight = dTOrig.LocalHeight();
        dTReplacement.AlignWith( dTOrig );
        dTReplacement.ResizeTo( dTOrig.Height(), 1 );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            dTReplacement.SetLocalEntry( iLocal, 0, (F)1 );
        LT.SetDiagonal( dTReplacement );
    }

    // Get ready for the local multiply
    DistMatrix<F,MR,STAR> XT_MR_STAR(g);
    XT_MR_STAR.AlignWith( LT );
    XT_MR_STAR = XT;

    // ZT[MC,* ] := LT[MC,MR] XT[MR,* ], 
    // XT[VC,* ].SumScatterFrom( ZT[MC,* ] )
    {
        DistMatrix<F,MC,STAR> ZT_MC_STAR(g);
        ZT_MC_STAR.AlignWith( LT );
        ZT_MC_STAR.ResizeTo( LT.Height(), XT.Width() );
        elem::internal::LocalGemm
        ( NORMAL, NORMAL, (F)1, LT, XT_MR_STAR, (F)0, ZT_MC_STAR );
        XT.SumScatterFrom( ZT_MC_STAR );
    }

    if( diag == UNIT )
        LT.SetDiagonal( dTOrig );

    if( LB.Height() != 0 )
    {
        // Set up for the local multiply
        XT_MR_STAR = XT;

        // ZB[MC,* ] := LB[MC,MR] XT[MR,* ]
        DistMatrix<F,MC,STAR> ZB_MC_STAR(g);
        ZB_MC_STAR.AlignWith( LB );
        ZB_MC_STAR.ResizeTo( LB.Height(), XT.Width() );
        elem::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, LB, XT_MR_STAR, (F)0, ZB_MC_STAR );
        XB.SumScatterUpdate( (F)1, ZB_MC_STAR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag, 
  DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontFastLowerBackwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
    if( orientation == NORMAL )
        throw std::logic_error("This solve must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();
    const int commSize = g.Size();
    const int commRank = g.VCRank();
    if( commSize == 1 )
    {
        LocalFrontLowerBackwardSolve
        ( orientation, diag, L.LockedLocalMatrix(), X.LocalMatrix() );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    elem::PartitionDown
    ( L, LT,
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    elem::PartitionDown
    ( X, XT,
         XB, snSize );

    // XT := XT - LB^{T/H} XB
    DistMatrix<F,STAR,STAR> Z( snSize, XT.Width(), g );
    if( XB.Height() != 0 )
    {
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, LB, XB, (F)0, Z );
        XT.SumScatterUpdate( (F)1, Z );
    }

    const int localTopHeight = LT.LocalHeight();
    std::vector<F> localDiag;
    if( diag == UNIT )
    {
        // Extract the diagonal of the top triangle and replace it with ones
        localDiag.resize( localTopHeight );
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            localDiag[iLocal] = LTBuffer[iLocal+i*LTLDim];
            LTBuffer[iLocal+i*LTLDim] = 1;
        }
    }

    // XT := LT^{T/H} XT
    elem::internal::LocalGemm( orientation, NORMAL, (F)1, LT, XT, (F)0, Z );
    XT.SumScatterFrom( Z );

    if( diag == UNIT )
    {
        // Put the diagonal back
        F* LTBuffer = LT.LocalBuffer();
        const int LTLDim = LT.LocalLDim();
        for( int iLocal=0; iLocal<localTopHeight; ++iLocal )
        {
            const int i = commRank + iLocal*commSize;
            LTBuffer[iLocal+i*LTLDim] = localDiag[iLocal];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag, 
  DistMatrix<F,MC,MR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontFastLowerBackwardSolve");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("This solve must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();
    const int commSize = g.Size();
    const int commRank = g.VCRank();
    if( commSize == 1 )
    {
        LocalFrontLowerBackwardSolve
        ( orientation, diag, L.LockedLocalMatrix(), X.LocalMatrix() );
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F,MC,MR> LT(g),
                        LB(g);
    elem::PartitionDown
    ( L, LT,
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    elem::PartitionDown
    ( X, XT,
         XB, snSize );

    // ZT[MR,* ] := -(LB[MC,MR])^{T/H} XB[MC,* ]
    // ZT[VR,* ].SumScatterFrom( ZT[MR,* ] )
    // ZT[VC,* ] := ZT[VR,* ]
    // XT[VC,* ] += ZT[VC,* ]
    DistMatrix<F,MR,STAR> ZT_MR_STAR( g );
    DistMatrix<F,VR,STAR> ZT_VR_STAR( g );
    ZT_MR_STAR.AlignWith( LB );
    ZT_MR_STAR.ResizeTo( snSize, XT.Width() );
    if( XB.Height() != 0 )
    {
        DistMatrix<F,MC,STAR> XB_MC_STAR( g );
        XB_MC_STAR.AlignWith( LB );
        XB_MC_STAR = XB;
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, LB, XB_MC_STAR, (F)0, ZT_MR_STAR );

        ZT_VR_STAR.SumScatterFrom( ZT_MR_STAR );
        DistMatrix<F,VC,STAR> ZT_VC_STAR( g );
        ZT_VC_STAR.AlignWith( XT );
        ZT_VC_STAR = ZT_VR_STAR;
        elem::Axpy( (F)1, ZT_VC_STAR, XT );
    }

    DistMatrix<F,MD,STAR> dTOrig(g), dTReplacement(g);
    if( diag == UNIT )
    {
        // Extract the diagonal of the top triangle and replace it with ones
        LT.GetDiagonal( dTOrig );
        dTReplacement.AlignWith( dTOrig );
        dTReplacement.ResizeTo( dTOrig.Height(), 1 );
        const int localHeight = dTOrig.LocalHeight();
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            dTReplacement.SetLocalEntry( iLocal, 0, (F)1 );
        LT.SetDiagonal( dTReplacement );
    }

    // ZT[MR,* ] := (LT[MC,MR])^{T/H} XT[MC,* ]
    // ZT[VR,* ].SumScatterFrom( ZT[MR,* ] )
    // XT[VC,* ] := ZT[VR,* ]
    {
        DistMatrix<F,MC,STAR> XT_MC_STAR( g );
        XT_MC_STAR.AlignWith( LT );
        XT_MC_STAR = XT;
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)1, LT, XT_MC_STAR, (F)0, ZT_MR_STAR );
        ZT_VR_STAR.SumScatterFrom( ZT_MR_STAR );
        XT = ZT_VR_STAR;
    }

    if( diag == UNIT )
        LT.SetDiagonal( dTOrig );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq

template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag, 
  DistMatrix<float,VC,STAR>& L,
  DistMatrix<float,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag, 
  DistMatrix<float,MC,MR  >& L,
  DistMatrix<float,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<float,VC,STAR>& L,
  DistMatrix<float,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<float,MC,MR  >& L,
  DistMatrix<float,VC,STAR>& X );

template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<double,VC,STAR>& L, 
  DistMatrix<double,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<double,MC,MR  >& L, 
  DistMatrix<double,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<double,VC,STAR>& L,
  DistMatrix<double,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<double,MC,MR  >& L,
  DistMatrix<double,VC,STAR>& X );

template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<Complex<float>,VC,STAR>& L, 
  DistMatrix<Complex<float>,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<Complex<float>,MC,MR  >& L, 
  DistMatrix<Complex<float>,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<Complex<float>,VC,STAR>& L, 
  DistMatrix<Complex<float>,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<Complex<float>,MC,MR  >& L, 
  DistMatrix<Complex<float>,VC,STAR>& X );

template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<Complex<double>,VC,STAR>& L, 
  DistMatrix<Complex<double>,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerForwardSolve
( Diagonal diag,
  DistMatrix<Complex<double>,MC,MR  >& L, 
  DistMatrix<Complex<double>,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<Complex<double>,VC,STAR>& L,
  DistMatrix<Complex<double>,VC,STAR>& X );
template void cliq::numeric::DistFrontFastLowerBackwardSolve
( Orientation orientation, Diagonal diag,
  DistMatrix<Complex<double>,MC,MR  >& L,
  DistMatrix<Complex<double>,VC,STAR>& X );
