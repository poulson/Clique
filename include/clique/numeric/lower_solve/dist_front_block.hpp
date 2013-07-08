/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename F>
void FrontBlockLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void FrontBlockLowerForwardSolve
( const DistMatrix<F>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void FrontBlockLowerBackwardSolve
( Orientation orientation, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X );

template<typename F>
void FrontBlockLowerBackwardSolve
( Orientation orientation, const DistMatrix<F>& L, DistMatrix<F,VC,STAR>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename F>
inline void FrontBlockLowerForwardSolve
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontBlockLowerForwardSolve");
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
    if( g.Size() == 1 )
    {
        FrontBlockLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    LockedPartitionDown
    ( L, LT, 
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    // XT := inv(ATL) XT
    DistMatrix<F,STAR,STAR> XT_STAR_STAR( XT );
    elem::LocalGemm( NORMAL, NORMAL, F(1), LT, XT_STAR_STAR, F(0), XT );

    // XB := XB - LB XT
    if( LB.Height() != 0 )
    {
        XT_STAR_STAR = XT;
        elem::LocalGemm( NORMAL, NORMAL, F(-1), LB, XT_STAR_STAR, F(1), XB );
    }
}

template<typename F>
inline void FrontBlockLowerForwardSolve
( const DistMatrix<F>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontBlockLowerForwardSolve");
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
    if( g.Size() == 1 )
    {
        FrontBlockLowerForwardSolve( L.LockedMatrix(), X.Matrix() );
        return;
    }

    // Separate the top and bottom portions of X and L
    const int snSize = L.Width();
    DistMatrix<F> LT(g),
                  LB(g);
    LockedPartitionDown
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

    {
        // ZT[MC,* ] := inv(ATL)[MC,MR] XT[MR,* ], 
        XT_MR_STAR = XT;
        DistMatrix<F,MC,STAR> ZT_MC_STAR(g);
        ZT_MC_STAR.AlignWith( LT );
        elem::LocalGemm( NORMAL, NORMAL, F(1), LT, XT_MR_STAR, ZT_MC_STAR );

        // XT[VC,* ] := SumScatterFrom( ZT[MC,* ] )
        XT.SumScatterFrom( ZT_MC_STAR );
    }

    if( LB.Height() != 0 )
    {
        // ZB[MC,* ] := LB[MC,MR] XT[MR,* ]
        XT_MR_STAR = XT;
        DistMatrix<F,MC,STAR> ZB_MC_STAR(g);
        ZB_MC_STAR.AlignWith( LB );
        elem::LocalGemm( NORMAL, NORMAL, F(1), LB, XT_MR_STAR, ZB_MC_STAR );

        // XB[VC,* ] -= ZB[MC,* ] = LB[MC,MR] XT[MR,* ]
        XB.SumScatterUpdate( F(-1), ZB_MC_STAR );
    }
}

template<typename F>
inline void FrontBlockLowerBackwardSolve
( Orientation orientation, 
  const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontBlockLowerBackwardSolve");
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
    if( g.Size() == 1 )
    {
        FrontBlockLowerBackwardSolve
        ( orientation, L.LockedMatrix(), X.Matrix() );
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F,VC,STAR> LT(g),
                          LB(g);
    LockedPartitionDown
    ( L, LT,
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    if( XB.Height() == 0 )
        return;

    // YT := LB^{T/H} XB
    DistMatrix<F,STAR,STAR> Z( g );
    elem::LocalGemm( orientation, NORMAL, F(1), LB, XB, Z );
    DistMatrix<F,VC,STAR> YT(g);
    YT.AlignWith( XT );
    YT.SumScatterFrom( Z );

    // XT := XT - inv(ATL) YT
    elem::LocalGemm( NORMAL, NORMAL, F(1), LT, YT, Z );
    XT.SumScatterUpdate( F(-1), Z );
}

template<typename F>
inline void FrontBlockLowerBackwardSolve
( Orientation orientation, const DistMatrix<F>& L, DistMatrix<F,VC,STAR>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontBlockLowerBackwardSolve");
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
    if( g.Size() == 1 )
    {
        FrontBlockLowerBackwardSolve
        ( orientation, L.LockedMatrix(), X.Matrix() );
        return;
    }

    const int snSize = L.Width();
    DistMatrix<F> LT(g),
                  LB(g);
    LockedPartitionDown
    ( L, LT,
         LB, snSize );
    DistMatrix<F,VC,STAR> XT(g),
                          XB(g);
    PartitionDown
    ( X, XT,
         XB, snSize );

    if( XB.Height() == 0 )
        return;


    DistMatrix<F,MR,STAR> ZT_MR_STAR( g );
    DistMatrix<F,VR,STAR> ZT_VR_STAR( g );
    ZT_MR_STAR.AlignWith( LB );
    DistMatrix<F,VC,STAR> YT(g);
    YT.AlignWith( XT );
    {
        // ZT[MR,* ] := (LB[MC,MR])^{T/H} XB[MC,* ]
        DistMatrix<F,MC,STAR> XB_MC_STAR( g );
        XB_MC_STAR.AlignWith( LB );
        XB_MC_STAR = XB;
        elem::LocalGemm
        ( orientation, NORMAL, F(1), LB, XB_MC_STAR, ZT_MR_STAR );

        // ZT[VR,* ].SumScatterFrom( ZT[MR,* ] )
        ZT_VR_STAR.SumScatterFrom( ZT_MR_STAR );

        // YT[VC,* ] := ZT[VR,* ]
        YT = ZT_VR_STAR;
    }

    {
        // ZT[MR,* ] := inv(ATL)[MC,MR] YT[MC,* ]
        DistMatrix<F,MC,STAR> YT_MC_STAR( g );
        YT_MC_STAR.AlignWith( LT );
        YT_MC_STAR = YT;
        elem::LocalGemm
        ( orientation, NORMAL, F(1), LT, YT_MC_STAR, ZT_MR_STAR );

        // ZT[VR,* ].SumScatterFrom( ZT[MR,* ] )
        ZT_VR_STAR.SumScatterFrom( ZT_MR_STAR );

        // ZT[VC,* ] := ZT[VR,* ]
        DistMatrix<F,VC,STAR> ZT_VC_STAR( g );
        ZT_VC_STAR.AlignWith( XT );
        ZT_VC_STAR = ZT_VR_STAR;

        // XT[VC,* ] -= ZT[VC,* ]
        elem::Axpy( F(-1), ZT_VC_STAR, XT );
    }
}

} // namespace cliq