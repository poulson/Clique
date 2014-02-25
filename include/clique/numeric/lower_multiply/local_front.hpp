/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
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
( int diagOff, const Matrix<T>& L, Matrix<T>& X, bool conjugate=false );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename T> 
void ModifyForTrmm
( Matrix<T>& D, int diagOff, std::vector<Matrix<T>>& diagonals )
{
    DEBUG_ONLY(CallStackEntry cse("ModifyForTrmm"))
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
  const std::vector<Matrix<T>>& diagonals )
{
    DEBUG_ONLY(CallStackEntry cse("ReplaceAfterTrmm"))
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
    DEBUG_ONLY(CallStackEntry cse("FrontLowerMultiply"))
    if( orientation == NORMAL )
        FrontLowerMultiplyNormal( diagOff, L, X );
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        FrontLowerMultiplyTranspose( diagOff, L, X, conjugate );
    }
}

template<typename T>
inline void FrontLowerMultiplyNormal
( int diagOff, const Matrix<T>& L, Matrix<T>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerMultiplyNormal");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal multiply:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( diagOff > 0 )
            LogicError("Diagonal offsets cannot be positive");
    )
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
        std::vector<Matrix<T>> diagonals;
        internal::ModifyForTrmm( LT, diagOff, diagonals );
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOff, diagonals );
    }
}

template<typename T>
inline void FrontLowerMultiplyTranspose
( int diagOff, const Matrix<T>& L, Matrix<T>& X, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FrontLowerMultiplyTranspose");
        if( L.Height() < L.Width() || L.Height() != X.Height() )
            LogicError
            ("Nonconformal solve:\n",
             DimsString(L,"L"),"\n",DimsString(X,"X"));
        if( diagOff > 0 )
            LogicError("Diagonal offsets cannot be positive");
    )
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
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
        std::vector<Matrix<T>> diagonals;
        internal::ModifyForTrmm( LT, diagOff, diagonals );
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOff, diagonals );
    }

    elem::Gemm( orientation, NORMAL, T(1), LB, XB, T(1), XT );
}

} // namespace cliq

#endif // ifndef CLIQ_NUMERIC_LOWERMULTIPLY_LOCALFRONT_HPP
