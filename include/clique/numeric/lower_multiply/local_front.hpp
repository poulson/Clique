/*
   Copyright (c) 2009-2013, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, and Stanford University
   All rights reserved.
 
   This file is part of Clique and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace cliq {

template<typename T>
void FrontLowerMultiply
( Orientation orientation, int diagOffset, const Matrix<T>& L, Matrix<T>& X );

template<typename T>
void FrontLowerMultiplyNormal
( int diagOffset, const Matrix<T>& L, Matrix<T>& X );

template<typename T>
void FrontLowerMultiplyTranspose
( Orientation orientation, int diagOffset, const Matrix<T>& L, Matrix<T>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename T> 
void ModifyForTrmm
( Matrix<T>& D, int diagOffset, std::vector<Matrix<T> >& diagonals )
{
#ifndef RELEASE
    CallStackEntry entry("ModifyForTrmm");
#endif
    diagonals.resize( -diagOffset );
    for( int i=0; i<-diagOffset; ++i )
        diagonals[i].ResizeTo( D.DiagonalLength(-i), 1 );

    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        for( int i=0; i<length; ++i )    
        {
            diagonals[i].Set( j, 0, D.Get(j+i,j) );
            D.Set( j+i, j, T(0) );
        }
    }
}

template<typename T> 
void ReplaceAfterTrmm
( Matrix<T>& D, int diagOffset, 
  const std::vector<Matrix<T> >& diagonals )
{
#ifndef RELEASE
    CallStackEntry entry("ReplaceAfterTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        for( int i=0; i<length; ++i )    
            D.Set( j+i, j, diagonals[i].Get(j,0) );
    }
}

} // namespace internal

template<typename T>
inline void FrontLowerMultiply
( Orientation orientation, int diagOffset,
  const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        FrontLowerMultiplyNormal( diagOffset, L, X );
    else
        FrontLowerMultiplyTranspose( orientation, diagOffset, L, X );
}

template<typename T>
inline void FrontLowerMultiplyNormal
( int diagOffset, const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontLowerMultiplyNormal");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal multiply:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
#endif
    // Danger, Will Robinson!
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT,
              LB;
    PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<T> XT, 
              XB;
    PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Gemm( NORMAL, NORMAL, T(1), LB, XT, T(1), XB );

    if( diagOffset == 0 )
    {
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T> > diagonals;
        internal::ModifyForTrmm( LT, diagOffset, diagonals );
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOffset, diagonals );
    }
}

template<typename T>
inline void FrontLowerMultiplyTranspose
( Orientation orientation, int diagOffset, const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("FrontLowerMultiplyTranspose");
    if( L.Height() < L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal solve:\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transposed");
    if( diagOffset > 0 )
        throw std::logic_error("Diagonal offsets cannot be positive");
#endif
    // Danger, Will Robinson!
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT,
              LB;
    PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<T> XT, 
              XB;
    PartitionDown
    ( X, XT,
         XB, L.Width() );

    if( diagOffset == 0 )
    {
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T> > diagonals;
        internal::ModifyForTrmm( LT, diagOffset, diagonals );
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diagOffset, diagonals );
    }

    elem::Gemm( orientation, NORMAL, T(1), LB, XB, T(1), XT );
}

} // namespace cliq
