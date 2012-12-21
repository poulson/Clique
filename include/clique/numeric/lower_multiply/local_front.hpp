/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T>
void LocalFrontLowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<T>& L, Matrix<T>& X );

template<typename T>
void LocalFrontLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset, const Matrix<T>& L, Matrix<T>& X );

template<typename T>
void LocalFrontLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<T>& L, Matrix<T>& X );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace internal {
using namespace elem;

template<typename T> 
void ModifyForTrmm
( Matrix<T>& D, UnitOrNonUnit diag, int diagOffset, 
  std::vector<Matrix<T> >& diagonals )
{
#ifndef RELEASE
    PushCallStack("ModifyForTrmm");
#endif
    if( diag == UNIT )
    {
        diagonals.resize( 1-diagOffset );
        for( int i=0; i<-diagOffset; ++i )
            diagonals[i].ResizeTo( D.DiagonalLength(-i), 1 );
        diagonals[-diagOffset].ResizeTo( D.DiagonalLength(-diagOffset), 1 );
    }
    else
    {
        diagonals.resize( -diagOffset );
        for( int i=0; i<-diagOffset; ++i )
            diagonals[i].ResizeTo( D.DiagonalLength(-i), 1 );
    }

    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        for( int i=0; i<length; ++i )    
        {
            diagonals[i].Set( j, 0, D.Get(j+i,j) );
            D.Set( j+i, j, T(0) );
        }
        if( diag == UNIT && j-diagOffset < height )
        {
            diagonals[-diagOffset].Set( j, 0, D.Get(j-diagOffset,j) );
            D.Set( j-diagOffset, j, T(1) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
void ReplaceAfterTrmm
( Matrix<T>& D, UnitOrNonUnit diag, int diagOffset, 
  const std::vector<Matrix<T> >& diagonals )
{
#ifndef RELEASE
    PushCallStack("ReplaceAfterTrmm");
#endif
    const int height = D.Height();
    for( int j=0; j<height; ++j )
    {
        const int length = std::min(-diagOffset,height-j);
        for( int i=0; i<length; ++i )    
            D.Set( j+i, j, diagonals[i].Get(j,0) );
        if( diag == UNIT && j-diagOffset < height )
            D.Set( j-diagOffset, j, diagonals[-diagOffset].Get(j,0) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void LocalFrontLowerMultiply
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("LocalFrontLowerMultiply");
#endif
    if( orientation == NORMAL )
        LocalFrontLowerMultiplyNormal( diag, diagOffset, L, X );
    else
        LocalFrontLowerMultiplyTranspose( orientation, diag, diagOffset, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void LocalFrontLowerMultiplyNormal
( UnitOrNonUnit diag, int diagOffset, const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("LocalFrontLowerMultiplyNormal");
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
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT,
              LB;
    elem::PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<T> XT, 
              XB;
    elem::PartitionDown
    ( X, XT,
         XB, L.Width() );

    elem::Gemm( NORMAL, NORMAL, T(1), LB, XT, T(1), XB );

    if( diagOffset == 0 )
    {
        elem::Trmm( LEFT, LOWER, NORMAL, diag, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T> > diagonals;
        internal::ModifyForTrmm( LT, diag, diagOffset, diagonals );
        elem::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diag, diagOffset, diagonals );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void LocalFrontLowerMultiplyTranspose
( Orientation orientation, UnitOrNonUnit diag, int diagOffset,
  const Matrix<T>& L, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("LocalFrontLowerMultiplyTranspose");
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
    Matrix<T>* LMod = const_cast<Matrix<T>*>(&L);
    Matrix<T> LT,
              LB;
    elem::PartitionDown
    ( *LMod, LT,
             LB, L.Width() );

    Matrix<T> XT, 
              XB;
    elem::PartitionDown
    ( X, XT,
         XB, L.Width() );

    if( diagOffset == 0 )
    {
        elem::Trmm( LEFT, LOWER, orientation, diag, T(1), LT, XT );
    }
    else
    {
        std::vector<Matrix<T> > diagonals;
        internal::ModifyForTrmm( LT, diag, diagOffset, diagonals );
        elem::Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), LT, XT );
        internal::ReplaceAfterTrmm( LT, diag, diagOffset, diagonals );
    }

    elem::Gemm( orientation, NORMAL, T(1), LB, XB, T(1), XT );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cliq
