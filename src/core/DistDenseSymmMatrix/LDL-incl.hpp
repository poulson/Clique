/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2010-2011 Jack Poulson <jack.poulson@gmail.com>
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

template<typename F>
void
clique::DistDenseSymmMatrix<F>::LocalLDL
( bool conjugate, int n, F* A, int lda )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LocalLDL");
    if( lda < n )
        throw std::logic_error
        ("Leading dimension cannot be smaller than height");
#endif
    if( n <= 1 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    std::vector<F> s21(n-1);
    for( int j=0; j<n; ++j )
    {
        const int a21Height = n - (j+1);

        // Extract the diagonal entry
        const F alpha11 = A[j+j*lda];    

        // Make a copy of a21 in s21 before scaling
        std::memcpy( &s21[0], &A[(j+1)+j*lda], a21Height*sizeof(F) );

        // a21 := a21 / alpha11
        const F alpha11Inv = static_cast<F>(1)/alpha11;
        {
            F* RESTRICT a21 = &A[(j+1)+j*lda];
            for( int i=0; i<a21Height; ++i )
                a21[i] *= alpha11Inv;
        }

        // A22 := A22 - s21 a21^[T/H]
        if( conjugate )
        {
            const F* RESTRICT a21 = &A[(j+1)+j*lda];
            for( int k=0; k<a21Height; ++k )
            {
                const F conjAlpha = Conj(a21[k]);
                F* RESTRICT A22Col = &A[(j+1)+(j+1+k)*lda];
                for( int i=k; i<a21Height; ++i )
                    A22Col[i] -= s21[i]*conjAlpha;
            }
        }
        else
        {
            const F* RESTRICT a21 = &A[(j+1)+j*lda];
            for( int k=0; k<a21Height; ++k )
            {
                const F alpha = a21[k];
                F* RESTRICT A22Col = &A[(j+1)+(j+1+k)*lda];
                for( int i=k; i<a21Height; ++i )
                    A22Col[i] -= s21[i]*alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::LDL( bool conjugate )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LDL");
#endif
    // We will frequently use these parameters, so let's use shorter names
    const int r = gridHeight_;
    const int c = gridWidth_;
    const int p = r*c;
    const int n = height_;
    const int nb = blockSize_;

    const int remainder = n % nb;
    const int numBlockCols = ( remainder==0 ? n/nb : n/nb+1 );

    const int maxA11Size = nb*nb;
    const int maxPackedA11Size = (nb*nb+nb)/2;
    std::vector<F> A11(maxA11Size), packedA11(maxPackedA11Size);

    const int maxA21LocalHeight = n/(r*nb) + nb;
    std::vector<F> diagAndA21( (maxA21LocalHeight+1)*nb );

    const int maxA21_VC_STAR_LocalHeight = n/(p*nb) + nb;
    std::vector<F> A21_VC_STAR( maxA21_VC_STAR_LocalHeight*nb );

    int jLocalBlock = 0;
    for( int jBlock=0; jBlock<numBlockCols; ++jBlock )
    {
        const int j = jBlock*nb;
        const int b = std::min(nb,n-j);
        const int ownerRow = jBlock % r;
        const int ownerCol = jBlock % c;
        const bool myRow = (ownerRow == gridRow_);
        const bool myCol = (ownerCol == gridCol_);

        const int iLocalBlock = (jBlock-gridRow_+r-1)/r;
        const int blockColHeight = localHeight_ - iLocalBlock*nb;
        const int A21LocalHeight = 
            ( myRow ? blockColHeight-b : blockColHeight );

        if( myCol )
        {
            F* blockCol = blockColBuffers_[jLocalBlock];
            F* A21 = ( myRow ? &blockCol[b] : blockCol );
            const int blockColLDim = blockColHeight;

            if( myRow )
            {
                // Perform the in-place LDL^[T/H] factorization
                LocalLDL( conjugate, b, blockCol, blockColLDim );

                // Pack the lower-triangle
                for( int k=0,offset=0; k<b; offset+=b-k,++k )
                    std::memcpy
                    ( &packedA11[offset], &blockCol[k+k*blockColLDim], 
                      (b-k)*sizeof(F) );
            }

            // Broadcast A11[MC,MR] from the owning row to form A11[* ,MR]
            const int packedSize = (b*b+b)/2;
            mpi::Broadcast( &packedA11[0], packedSize, ownerRow, colComm_ );

            // Unpack A11
            for( int k=0,offset=0; k<b; offset+=b-k,++k )
                std::memcpy( &A11[k+k*b], &packedA11[offset], (b-k)*sizeof(F) );

            // A21 := A21 trilu(A11)^-[T/H] 
            const char option = ( conjugate ? 'C' : 'T' );
            blas::Trsm
            ( 'R', 'L', option, 'U', A21LocalHeight, b, 
              (F)1, &A11[0], b, A21, blockColLDim );

            // Copy D11[* ,MR] and A21[MC,MR] into a buffer
            for( int j=0; j<b; ++j )
                diagAndA21[j] = A11[j+j*b];
            for( int j=0; j<b; ++j )
                std::memcpy
                ( &diagAndA21[b+j*A21LocalHeight], &A21[j*blockColLDim], 
                  A21LocalHeight*sizeof(F) );
        }

        // Broadcast D11[* ,MR] and A21[MC,MR] within rows to form 
        // D11[* ,* ] and A21[MC,* ]
        mpi::Broadcast
        ( &diagAndA21[0], (A21LocalHeight+1)*b, ownerCol, rowComm_ );

        if( gridHeight_ == gridWidth_ )
        {
            if( gridCol_ == gridRow_ )
            {
                // Just perform a memcpy
                // TODO
            }
            else
            {
                // Pairwise Send/Recv to form S21[MR,* ] from A21[MC,* ]
                // TODO
            }
        }
        else
        {
            // Locally copy A21[VC,* ] out of A21[MC,* ]
            const int A21_VC_STAR_Align = (jBlock+1) % p;
            //const int myFirstA21Block = ?

            // TODO

            // SendRecv permutation to form A21[VR,* ] from A21[VC,* ]
            // TODO

            // Gather A21[VR,* ] within rows to form S21[MR,* ]
            // TODO
        }

        // S21[MC,* ] := A21[MC,* ] D11^-1
        // TODO

        // A22[MC,MR] -= S21[MC,* ] A21[MR,* ]^[T/H]
        // TODO

        if( myCol )
            ++jLocalBlock;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::LDLT()
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LDLT");
#endif
    LDL( false );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::LDLH()
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LDLH");
#endif
    LDL( true );
#ifndef RELEASE
    PopCallStack();
#endif
}

