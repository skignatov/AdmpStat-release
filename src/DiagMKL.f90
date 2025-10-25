Subroutine DiagMKL(n,A,lda,W,V,info,isort)
    Use lapack95
    Implicit Real(8) (A-H,O-Z)

!    Include 'mkl.fi'  		! Can be used instead of "USE lapack95"
    
    Real(8) A(lda,n),W(n),Work(3*n-1),V(lda,n)
    Character(1) jobz,uplo

! Program for real symmetric matrix diagonalization using the MKL (really LAPACK) library routine dsyev (see: http://cali2.unilim.fr/intel-xe/mkl/mklman)
! solving the simple diagonalization problem Ax=ex (for general diag problem Ax=eBx see the routine dsygv)
!
! (Please take into account thet MKL contains also modern diagonalization routine dfeast_syev and its modifications 
! solving also generalized diagonalization problem Ax=eBx for any types of matrices by innovative robust method FEAST)
! 
!
! Input:
!   n - problem size
!   A(lda,n)    - matrix to be diagonalized
!   lda         - leading dimension
!   isort       - Sort method for eigenvalues: if isort<0 - descending order, otherwise -ascending order
!
! Output:
!   W(1:n)      - Eigenvalues (normally sorted in ascending order by ?syev)
!   V(1:n,1:n)  - Eigenvectors (by columns)
!   info        - if 0 - normal termination, otherwise - error described in ?syev documentation
    
    
    V(1:n,1:n)=A(1:n,1:n)
    
    jobz='V'
    uplo='U'
    lwork=3*n-1
    call dsyev(jobz, uplo, n, V, lda, w, work, lwork, info)     ! Diagonalization. Eigenvalues are in ascending order

    If (isort<0) Then   ! Sorting eigenvalues in descending order
        n1=n+1
        Do i=1,n
            k=n1-i
            Work(k)=W(i)
        Enddo
        W(1:n)=Work(1:n)
        Do j=1,n
            Do i=1,n
                k=n1-i
                Work(k)=V(j,i)
            Enddo
            V(j,1:n)=Work(1:n)
        Enddo
    Endif
    
End
