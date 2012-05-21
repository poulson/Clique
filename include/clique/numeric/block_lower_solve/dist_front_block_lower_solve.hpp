/*
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
#ifndef CLIQUE_NUMERIC_DIST_FRONT_BLOCK_LOWER_SOLVE_HPP
#define CLIQUE_NUMERIC_DIST_FRONT_BLOCK_LOWER_SOLVE_HPP 1

namespace cliq {
namespace numeric {

template<typename F>
void DistFrontBlockLowerForwardSolve
( DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void DistFrontBlockLowerForwardSolve
( DistMatrix<F,MC,MR>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void DistFrontBlockLowerBackwardSolve
( Orientation orientation, DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void DistFrontBlockLowerBackwardSolve
( Orientation orientation, DistMatrix<F,MC,MR>& L, DistMatrix<F,VC,STAR>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void DistFrontLowerForwardSolve
( DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontBlockLowerForwardSolve");
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
        numeric::LocalFrontBlockLowerForwardSolve
        ( L.LockedLocalMatrix(), X.LocalMatrix() );
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

    // Get a copy of all of XT
    DistMatrix<F,STAR,STAR> XT_STAR_STAR( XT );

    // XT := LT XT
    elem::internal::LocalGemm
    ( NORMAL, NORMAL, (F)1, LT, XT_STAR_STAR, (F)0, XT );

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
inline void DistFrontBlockLowerForwardSolve
( DistMatrix<F,MC,MR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontBlockLowerForwardSolve");
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
        numeric::LocalFrontBlockLowerForwardSolve
        ( L.LockedLocalMatrix(), X.LocalMatrix() );
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

    // Get ready for the local multiply
    DistMatrix<F,MR,STAR> XT_MR_STAR(g);
    XT_MR_STAR.AlignWith( LT );
    XT_MR_STAR = XT;

    // ZT[MC,* ] := inv(ATL)[MC,MR] XT[MR,* ], 
    // XT[VC,* ].SumScatterFrom( ZT[MC,* ] )
    {
        DistMatrix<F,MC,STAR> ZT_MC_STAR(g);
        ZT_MC_STAR.AlignWith( LT );
        ZT_MC_STAR.ResizeTo( LT.Height(), XT.Width() );
        elem::internal::LocalGemm
        ( NORMAL, NORMAL, (F)1, LT, XT_MR_STAR, (F)0, ZT_MC_STAR );
        XT.SumScatterFrom( ZT_MC_STAR );
    }

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
inline void DistFrontBlockLowerBackwardSolve
( Orientation orientation, UnitOrNonUnit diag, 
  DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontBlockLowerBackwardSolve");
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
        LocalFrontBlockLowerBackwardSolve
        ( orientation, L.LockedLocalMatrix(), X.LocalMatrix() );
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

    // YT := LB^{T/H} XB
    DistMatrix<F,VC,STAR> YT(g);
    YT.AlignWith( XT );
    Zeros( XT.Height(), XT.Width(), YT );
    DistMatrix<F,STAR,STAR> Z( snSize, XT.Width(), g );
    if( XB.Height() != 0 )
    {
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)1, LB, XB, (F)0, Z );
        YT.SumScatterFrom( Z );
    }

    // XT := XT - inv(ATL) YT
    elem::internal::LocalGemm( NORMAL, NORMAL, (F)1, LT, YT, (F)0, Z );
    XT.SumScatterUpdate( (F)-1, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void DistFrontBlockLowerBackwardSolve
( Orientation orientation, DistMatrix<F,MC,MR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    PushCallStack("numeric::DistFrontBlockLowerBackwardSolve");
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
        LocalFrontBlockLowerBackwardSolve
        ( orientation, L.LockedLocalMatrix(), X.LocalMatrix() );
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
    // YT[VC,* ] := ZT[VR,* ]
    DistMatrix<F,VC,STAR> YT(g);
    YT.AlignWith( XT );
    Zeros( XT.Height(), XT.Width(), YT );
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
        YT = ZT_VR_STAR;
    }

    // ZT[MR,* ] := -(LT[MC,MR])^{T/H} YT[MC,* ]
    // ZT[VR,* ].SumScatterFrom( ZT[MR,* ] )
    // ZT[VC,* ] := ZT[VR,* ]
    // XT[VC,* ] += ZT[VR,* ]
    {
        DistMatrix<F,MC,STAR> YT_MC_STAR( g );
        YT_MC_STAR.AlignWith( LT );
        YT_MC_STAR = YT;
        elem::internal::LocalGemm
        ( orientation, NORMAL, (F)-1, LT, YT_MC_STAR, (F)0, ZT_MR_STAR );
        ZT_VR_STAR.SumScatterFrom( ZT_MR_STAR );
        DistMatrix<F,VC,STAR> ZT_VC_STAR( g );
        ZT_VC_STAR.AlignWith( XT );
        ZT_VC_STAR = ZT_VR_STAR;
        elem::Axpy( (F)1, ZT_VC_STAR, XT );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace numeric
} // namespace cliq

#endif // CLIQUE_NUMERIC_DIST_FRONT_BLOCK_LOWER_SOLVE_HPP
