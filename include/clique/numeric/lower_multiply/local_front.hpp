/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CLIQ_NUMERIC_LOWERMULTIPLY_LOCALFRONT_HPP
#define CLIQ_NUMERIC_LOWERMULTIPLY_LOCALFRONT_HPP

namespace cliq {

template<typename T>
void FrontLowerMultiply
( Orientation orientation, int diagOff, const Matrix<T>& L, Matrix<T>& X );

template<typename T>
void FrontLowerMultiplyNormal( int diagOff, const Matrix<T>& L, Matrix<T>& X );

template<typename T>
void FrontLowerMultiplyTranspose
( Orientation orientation, int diagOff, const Matrix<T>& L, Matrix<T>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename T> 
void ModifyForTrmm
( Matrix<T>& D, int diagOff, std::vector<Matrix<T> >& diagonals )
{
#ifndef RELEASE
    CallStackEntry cse("ModifyForTrmm");
#endif
    diagonals.resize( -diagOff );
    for( int i=0; i<-diagOff; ++i )
        diagonals[i].ResizeTo( D.DiagonalLength(-i), 1 );

    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOff,height-j);
        for( int i=0; i<length; ++i )    
        {
            diagonals[i].Set( j, 0, D.Get(j+i,j) );
            D.Set( j+i, j, T(0) );
        }
    }
}

template<typename T> 
void ReplaceAfterTrmm
( Matrix<T>& D, int diagOff, 
  const std::vector<Matrix<T> >& diagonals )
{
#ifndef RELEASE
    CallStackEntry cse("ReplaceAfterTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOff,height-j);
        for( int i=0; i<length; ++i )    
            D.Set( j+i, j, diagonals[i].Get(j,0) );
    }
}

} // namespace internal

template<typename T>
inline void FrontLowerMultiply
( Orientation orientation, int diagOff, const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        FrontLowerMultiplyNormal( diagOff, L, X );
    else
        FrontLowerMultiplyTranspose( orientation, diagOff, L, X );
}

template<typename T>
inline void FrontLowerMultiplyNormal
( int diagOff, const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerMultiplyNormal");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( diagOff > 0 )
        LogicError("Diagonal offsets cannot be positive");
#endif
    // Danger, Will Robinson!
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT, LB;
    PartitionDown( *LMod, LT, LB, L.Width() );
    Matrix<T> XT, XB;
    PartitionDown( X, XT, XB, L.Width() );

    elem::Gemm( NORMAL, NORMAL, T(1), LB, XT, T(1), XB );

    if( diagOff == 0 )
    {
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T> > diagonals;
        internal::ModifyForTrmm( LT, diagOff, diagonals );
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOff, diagonals );
    }
}

template<typename T>
inline void FrontLowerMultiplyTranspose
( Orientation orientation, int diagOff, const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry cse("FrontLowerMultiplyTranspose");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        LogicError( msg.str() );
    }
    if( orientation == NORMAL )
        LogicError("Orientation must be (conjugate-)transposed");
    if( diagOff > 0 )
        LogicError("Diagonal offsets cannot be positive");
#endif
    // Danger, Will Robinson!
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT, LB;
    PartitionDown( *LMod, LT, LB, L.Width() );
    Matrix<T> XT, XB;
    PartitionDown( X, XT, XB, L.Width() );
    if( diagOff == 0 )
    {
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T> > diagonals;
        internal::ModifyForTrmm( LT, diagOff, diagonals );
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOff, diagonals );
    }

    elem::Gemm( orientation, NORMAL, T(1), LB, XB, T(1), XT );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERMULTIPLY_LOCALFRONT_HPP
