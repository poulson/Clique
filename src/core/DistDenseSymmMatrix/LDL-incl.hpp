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
clique::DistDenseSymmMatrix<F>::LDL( bool conjugate, int q )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LDL");
    if( q > height_ )
        throw std::logic_error("Short-circuit parameter was too large");
#endif
    const char option = ( conjugate ? 'C' : 'T' );

    const int r = gridHeight_;
    const int c = gridWidth_;
    const int p = r*c;
    const int MCRank = gridRow_;
    const int MRRank = gridCol_;
    const int VCRank = MCRank + r*MRRank;
    const int VRRank = MRRank + c*MCRank;

    const int n = height_;
    const int nb = blockSize_;
    const int numBlockCols = BlockLength( n, nb );
    const int mLocal = LocalLength( n, nb, 0, MCRank, r );
    const int nLocal = LocalLength( n, nb, 0, MRRank, c );

    const int maxA11Size = nb*nb;
    const int maxPackedA11Size = (nb*nb+nb)/2;
    std::vector<F> A11(maxA11Size), packedA11(maxPackedA11Size);

    const int maxA21_MC_STAR_LocalHeight = (n/(r*nb)+1)*nb;
    const int maxA21_MR_STAR_LocalHeight = (n/(c*nb)+1)*nb;
    std::vector<F> diagAndA21_MC_STAR( (maxA21_MC_STAR_LocalHeight+1)*nb );
    std::vector<F> A21_MR_STAR( maxA21_MR_STAR_LocalHeight*nb );

    const int maxA21VectorLocalHeight = (n/(p*nb)+1)*nb;
    std::vector<F> A21_VC_STAR, A21_VR_STAR, A21_MR_STAR_Temp;
    if( r != c )
    {
        A21_VC_STAR.resize( maxA21VectorLocalHeight*nb );
        A21_VR_STAR.resize( maxA21VectorLocalHeight*nb );
        A21_MR_STAR_Temp.resize( maxA21VectorLocalHeight*r*nb );
    }

    int jLocalBlock = 0;
    const int numFactorBlockCols = BlockLength( q, nb );
    for( int jBlock=0; jBlock<numFactorBlockCols; ++jBlock )
    {
        const int j = jBlock*nb;
        const int b = std::min(nb,n-j);
        const int bf = std::min(nb,q-j);
        const int bRem = b - bf;

        const int ownerRow = jBlock % r;
        const int ownerCol = jBlock % c;
        const bool myRow = (ownerRow == MCRank);
        const bool myCol = (ownerCol == MRRank);

        const int blockColLocalHeight = 
            LocalLength( n-j, nb, ownerRow, MCRank, r );
        const int A21_MC_STAR_LocalHeight = 
            LocalLength( n-j-b, nb, (jBlock+1)%r, MCRank, r );
        const int A21_MR_STAR_LocalHeight = 
            LocalLength( n-j-b, nb, (jBlock+1)%c, MRRank, c );
        const int A21_VC_STAR_LocalHeight = 
            LocalLength( n-j-b, nb, (jBlock+1)%p, VCRank, p );
        const int A21_VR_STAR_LocalHeight = 
            LocalLength( n-j-b, nb, (jBlock+1)%p, VRRank, p );

        // If we are in the owning process column, then we first factor/update 
        // it within our process column and then begin distributing it to other
        // process columns.
        if( myCol )
        {
            F* blockCol = blockColBuffers_[jLocalBlock];
            F* A21 = ( myRow ? &blockCol[b] : blockCol );
            const int blockColLDim = blockColLocalHeight;

            if( myRow )
            {
                // Perform the in-place LDL^[T/H] factorization on the top-left
                // bf x bf block
                LocalLDL( conjugate, bf, blockCol, blockColLDim );

                // Solve against the bottom portion of the diagonal block
                // Technically this is in A21.
                blas::Trsm
                ( 'R', 'L', option, 'U', bRem, bf,
                  (F)1, &blockCol[0], blockColLDim, 
                        &blockCol[bf], blockColLDim );

                // Pack the first bf columns of the lower triangle
                for( int t=0,offset=0; t<bf; offset+=b-t,++t )
                    std::memcpy
                    ( &packedA11[offset], &blockCol[t+t*blockColLDim], 
                      (b-t)*sizeof(F) );
           }

            // Broadcast A11[MC,MR] from the owning row to form A11[* ,MR]
            const int packedSize = (b*b+b-bRem*bRem-bRem)/2;
            mpi::Broadcast( &packedA11[0], packedSize, ownerRow, colComm_ );

            // Unpack A11
            for( int k=0,offset=0; k<bf; offset+=b-k,++k )
                std::memcpy
                ( &A11[k+k*b], &packedA11[offset], (b-k)*sizeof(F) );

            // Finish A21 := A21 trilu(A11)^-[T/H] 
            // (if bf != b, we already did this solve on the diag block)
            blas::Trsm
            ( 'R', 'L', option, 'U', A21_MC_STAR_LocalHeight, bf, 
              (F)1, &A11[0], b, A21, blockColLDim );

            // Copy D11[* ,MR] and A21[MC,MR] into a buffer
            for( int t=0; t<bf; ++t )
                diagAndA21_MC_STAR[t] = A11[t+t*b];
            for( int t=0; t<bf; ++t )
                std::memcpy
                ( &diagAndA21_MC_STAR[bf+t*A21_MC_STAR_LocalHeight], 
                  &A21[t*blockColLDim], A21_MC_STAR_LocalHeight*sizeof(F) );

            // Locally solve against D11 to finish forming the first bf columns
            // of this block column
            for( int t=0; t<bf; ++t )
            {
                const F invDelta = static_cast<F>(1)/A11[t+t*b];
                if( myRow )
                    for( int s=0; s<bRem; ++s )
                        blockCol[(bf+s)+t*blockColLDim] *= invDelta;
                F* A21Col = &A21[t*blockColLDim];
                for( int s=0; s<A21_MC_STAR_LocalHeight; ++s )
                    A21Col[s] *= invDelta;
            }

            // Update the last b-bf=bRem columns of the block column
            if( myRow )
            {
                // Update the diagonal block
                blas::Gemm
                ( 'N', option, bRem, bRem, bf,
                  (F)-1, &blockCol[bf], blockColLDim, &A11[bf], b,
                  (F)1,  &blockCol[bf+bf*blockColLDim], blockColLDim );
            }
            // Update the rest of the block column
            blas::Gemm
            ( 'N', option, A21_MC_STAR_LocalHeight, bRem, bf,
              (F)-1, A21, blockColLDim, &A11[bf], b,
              (F)1,  &A21[bf*blockColLDim], blockColLDim );

            ++jLocalBlock;
        }

        // Broadcast D11[* ,MR] and A21[MC,MR] within rows to form 
        // D11[* ,* ] and A21[MC,* ]
        mpi::Broadcast
        ( &diagAndA21_MC_STAR[0], (A21_MC_STAR_LocalHeight+1)*bf, 
          ownerCol, rowComm_ );
        F* A21_MC_STAR = &diagAndA21_MC_STAR[bf];

        if( r == c )
        {
            throw std::logic_error("Fast transpose not yet written");
            if( MCRank == MRRank )
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
            const int A21_VC_STAR_Shift = 
                BlockShift( n-j-b, b, (jBlock+1)%p, VCRank, p );
            const int A21_MC_STAR_Shift =
                BlockShift( n-j-b, b, (jBlock+1)%r, MCRank, r );
            const int A21_VC_STAR_RelativeShift = 
                (A21_VC_STAR_Shift-A21_MC_STAR_Shift)/r;
            const int A21_VC_STAR_LocalBlockHeight = 
                BlockLength( A21_VC_STAR_LocalHeight, b );
            for( int s=0; s<A21_VC_STAR_LocalBlockHeight; ++s )
            {
                const int iLocalBlock = A21_VC_STAR_RelativeShift + s*c;
                const F* origBlock = &A21_MC_STAR[iLocalBlock*b];
                F* copiedBlock = &A21_VC_STAR[s*b];
                const int blockHeight = 
                    std::min(b,A21_MC_STAR_LocalHeight-iLocalBlock*b);
                for( int t=0; t<bf; ++t )
                    std::memcpy
                    ( &copiedBlock[t*A21_VC_STAR_LocalHeight], 
                      &origBlock[t*A21_MC_STAR_LocalHeight], 
                      blockHeight*sizeof(F) );
            }

            // SendRecv permutation to form A21[VR,* ] from A21[VC,* ]
            const int sendVCRank = RelabelVRToVC( VCRank );
            const int recvVCRank = RelabelVCToVR( VCRank );
            const int A21_VR_STAR_LocalHeight = 
                LocalLength( n-(j+b), b, (jBlock+1)%p, VRRank, p );
            mpi::SendRecv
            ( &A21_VC_STAR[0], A21_VC_STAR_LocalHeight*bf, sendVCRank, 0,
              &A21_VR_STAR[0], A21_VR_STAR_LocalHeight*bf, recvVCRank, 0, 
              cartComm_ );

            // AllGather A21[VR,* ] within columns to form A21[MR,* ]
            const int portionSize = ((n-j-b)/(p*b)+1)*b*bf;
            mpi::AllGather
            ( &A21_VR_STAR[0],      portionSize, 
              &A21_MR_STAR_Temp[0], portionSize, colComm_ );
            const int A21_MR_STAR_Shift = 
                BlockShift( n-j-b, b, (jBlock+1)%c, MRRank, c );
            for( int k=0; k<r; ++k )
            {
                const int thisVRRank = MRRank + k*c;
                const int thisVRShift = 
                    BlockShift( n-j-b, b, (jBlock+1)%p, thisVRRank, p );
                const int thisVRLocalLength = 
                    LocalLength( n-j-b, b, (jBlock+1)%p, thisVRRank, p );
                const int thisVRLocalBlockLength = 
                    BlockLength( thisVRLocalLength, b );
                const int thisVRRelativeShift = 
                    (thisVRShift-A21_MR_STAR_Shift)/c;
                const F* thisData = &A21_MR_STAR_Temp[k*portionSize];

                for( int s=0; s<thisVRLocalBlockLength; ++s )
                {
                    const int blockHeight = std::min(b,thisVRLocalLength-s*b);
                    const F* tempBuf = &thisData[s*b];
                    F* finalBuf = &A21_MR_STAR[(thisVRRelativeShift+s*r)*b];
                    for( int t=0; t<bf; ++t )
                        std::memcpy
                        ( &finalBuf[t*A21_MR_STAR_LocalHeight], 
                          &tempBuf[t*thisVRLocalLength], 
                          blockHeight*sizeof(F) );
                }
            }
        }

        // S21[MC,* ] := A21[MC,* ] D11^-1
        for( int t=0; t<bf; ++t )
        {
            const F invDelta = static_cast<F>(1)/diagAndA21_MC_STAR[t];
            F* A21Col = &A21_MC_STAR[t*A21_MC_STAR_LocalHeight];
            for( int s=0; s<A21_MC_STAR_LocalHeight; ++s )
                A21Col[s] *= invDelta;
        }

        // A22[MC,MR] -= S21[MC,* ] A21[MR,* ]^[T/H]
        // TODO: Increase widths of gemms for better performance
        const int A21_MR_STAR_LocalBlockHeight = 
            BlockLength( A21_MR_STAR_LocalHeight, b );
        for( int t=0; t<A21_MR_STAR_LocalBlockHeight; ++t )
        {
            const int blockWidth = std::min(b,A21_MR_STAR_LocalHeight-t*b);
            const int blockColLocalHeight = 
                blockColLocalHeights_[jLocalBlock+t];

            const F* A21_MC_STAR_Sub = 
                &A21_MC_STAR[A21_MC_STAR_LocalHeight-blockColLocalHeight];
            const F* A21_MR_STAR_Sub = &A21_MR_STAR[t*b];
            F* blockCol = blockColBuffers_[jLocalBlock+t];

            blas::Gemm
            ( 'N', option, blockColLocalHeight, blockWidth, bf,
              (F)-1, A21_MC_STAR_Sub, A21_MC_STAR_LocalHeight,
                     A21_MR_STAR_Sub, A21_MR_STAR_LocalHeight,
              (F)1,  blockCol,        blockColLocalHeight );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::LDLT( int q )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LDLT");
#endif
    LDL( false, q );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
clique::DistDenseSymmMatrix<F>::LDLH( int q )
{
#ifndef RELEASE
    PushCallStack("DistDenseSymmMatrix::LDLH");
#endif
    LDL( true, q );
#ifndef RELEASE
    PopCallStack();
#endif
}

