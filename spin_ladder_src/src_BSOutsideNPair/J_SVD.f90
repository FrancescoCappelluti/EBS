!/************************************************************************
MODULE J_SVD
  USE PARAMETERS,   ONLY: dp,STDOUT,STDERR
  USE PRINT_TOOLS,  ONLY: write_matrix

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Calculate_J
!************************************************************************
  CONTAINS

    SUBROUTINE Calculate_J(energies,M_matrix,A,M_ave,J_values,N_state,N_Center,N_Pair,unit_out)
        REAL(kind=dp), INTENT(IN), DIMENSION(:)   :: energies
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:) :: M_matrix
        REAL(kind=dp), INTENT(OUT), DIMENSION(:)  :: M_ave
        REAL(kind=dp), INTENT(OUT), DIMENSION(:)  :: J_values
        INTEGER, INTENT(IN)                       :: N_state,N_Center,N_Pair
        INTEGER, INTENT(IN)                       :: unit_out
        REAL(kind=dp), INTENT(OUT), DIMENSION(N_state,N_Pair) :: A
        REAL(kind=dp), DIMENSION(N_Pair) :: S
        REAL(kind=dp), DIMENSION(N_state,N_Pair) :: U
        REAL(kind=dp), DIMENSION(N_Pair,N_Pair) :: Vt
        REAL(kind=dp), DIMENSION(N_Pair,N_Pair) :: Sigma,temp_ave

!        real(kind=dp),  :: A(:,:), U(:,:), Vt(:,:), Sigma(:,:)
!        real(kind=dp), ALLOCATABLE :: M_matrix(:,:),M_ave(:), Mul_Matrix(:)
!        real(kind=dp), ALLOCATABLE :: tmp(:),tmp_ava(:),temp_ave(:,:)
        real(kind=dp), ALLOCATABLE ::  work(:)
!        INTEGER, ALLOCATABLE :: desca(:),descu(:),descvt(:)
        INTEGER :: ifail, i,j,k,p, info, lwork
        real(kind=dp) :: fatt_au2cm1
!        INTEGER   ::    iostat
!        INTEGER, PARAMETER     :: unit_in=27,unit_out=29
!        CHARACTER(len=20)  :: FMT
        fatt_au2cm1=219474.63_dp

        write(unit_out,'(a)')'#     M MATRIX: M_(k,i) [k=BS State; i=Spin Center]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Center
        call write_matrix(matrix=TRANSPOSE(M_matrix),cutoff=N_Center,unit=unit_out)
        write(unit_out,*)

        ! Make averages along columns (absolute value)
        M_ave=0.0_dp
        do i=1,N_state
          do j=1,N_Center
            M_ave(j)=M_ave(j)+abs(M_matrix(i,j))
          enddo
        enddo
        M_ave(:)=M_ave(:)/N_state

        ! Make the proper multiplications to build the matrix A = M_i*M_j/4
        A=0.0_dp
        do k=1,N_state
          p=1
          do i=1,N_Center
            do j=i+1,N_Center
              A(k,p)=sign(M_ave(i),M_matrix(k,i))*sign(M_ave(j),M_matrix(k,j))/4.0_dp
              p=p+1
            enddo
          enddo
        enddo


        U=0.0_dp
        Vt=0.0_dp
        s=0.0_dp
        lwork=6*N_Pair
        allocate(work(lwork))

        lwork=(N_state)**2 ; ifail=0

        write(unit_out,'(a)')'#     A MATRIX: (M_i*M_j)_(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

!        do i=1,N_state
!          write (stdout,'(12f10.5)') (A(i,j),j=1,N_Pair)
!        enddo
        A=A*2
!        E=-(2N_state)*J

call dgesvd('S','S', N_state, N_Pair, A, N_state, s, U, N_state, Vt, N_Pair, work, lwork, info)


!Scrivo i risultati della procedura SVD
        write(unit_out,'(a)')'#     U MATRIX: U_(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(U),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

!        write (stdout,*) 'Matrix U:'
!        do i=1,N_state
!          write (stdout,'(12f10.5)') (U(i,j),j=1,N_Pair)
!        enddo
!        do i=1,N_Pair
!          write (stdout,'(12f10.5)') (Vt(i,j),j=1,N_Pair)
!        enddo

        write(unit_out,'(a)')'#    Vt MATRIX: Vt_(p1,p1) [p1,p2=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_Pair,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(Vt),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        write(unit_out,'(a)')'#    Sigma Diagonal Values:'
        write(unit_out,'(a,I5)')'# ',N_Pair
        write(unit_out,'(6F10.5)') (s(j),j=1,N_Pair)
        write(unit_out,*)

!        write (stdout,*) 'diag Matrix S:'
!        write (stdout,'(12f10.5)') (s(j),j=1,N_Pair)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Verifico moltiplicando tra di loro le matrici A=U*Sigma*Vt
!!! Verified by multiplying Matrix A=U*Sigma*Vt
!        Sigma=0.0_dp
!        do i=1,N_Pair
!          Sigma(i,i)=s(i)
!        enddo
!
!       temp_ave=0.0_dp
!       do p=1,N_Pair
!        do i=1,N_Pair
        !           do j=1,N_Pair
        
!            temp_ave(i,j)=temp_ave(i,j)+Sigma(i,p)*Vt(p,j)
!           enddo
!        enddo
!       enddo
!       A=0
!       do p=1,N_Pair
!        do i=1,N_state
!           do j=1,N_Pair
!            A(i,j)=A(i,j)+U(i,p)*temp_ave(p,j)
!           enddo
!        enddo
!       enddo
!       !Tolgo il raddoppio precedente
!       A=A/2
!       write (stdout,*) 'Matrix A: '
!       do i=1,N_state
!          write (stdout,'(12f10.5)') (A(i,j),j=1,N_Pair)
!       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! J CALCULATION
!        write (stdout,*) 'Starting J_ij calculations.'


       !Inverto la matrice diagonale
       !Matrix A=M*M, then calculate the inverse matrix A in following section
       !Inverse diagonalize terms
        Sigma=0.0_dp
        do i=1,N_Pair
          Sigma(i,i)=1/s(i)
        enddo

        !Calcolo la Matrice inversa di A: (A^-1)^T=U*Sigma^-1*Vt
        !get the inverse of Sigma Matrix(inverse the diagonalize terms)
        !then multiply by U and Vt in order to get the inverse matrix of A in following two paragraphs
       temp_ave=0.0_dp
       do p=1,N_Pair
        do i=1,N_Pair
           do j=1,N_Pair
            temp_ave(i,j)=temp_ave(i,j)+Sigma(i,p)*Vt(p,j)
           enddo
        enddo
       enddo
       A=0.0_dp
       do p=1,N_Pair
        do i=1,N_state
           do j=1,N_Pair
            A(i,j)=A(i,j)+U(i,p)*temp_ave(p,j)
           enddo
        enddo
        enddo


        J_values=0.0_dp
        do i=1,N_Pair
          do j=1,N_state
              J_values(i)=J_values(i)-A(j,i)*energies(j)
          enddo
        enddo

        !print J coupling values into svd_matrix.out
        write(unit_out,'(a)')'#    J coupling:'
        write(unit_out,'(a,I5)')'# ',N_Pair
        write(unit_out,'(ES14.7)') (J_values(j),j=1,N_Pair)
        write(unit_out,*)

        J_values=J_values*fatt_au2cm1

!        !Scrivo J
!        write (stdout,*) ' The J values: '
!        do i=1,N_Pair
!           write (stdout,'(1f15.5)') J_values(i)
!        enddo

!        open (unit=unit_out,file="SpinLadder.inp",action="write",status="replace")
!        write (unit_out,'(I3)')N_Center
!        write(FMT,'(a,i,a)')'(',N_Center,'F5.1)'
!        write (unit_out,FMT)(real(int(M_ave(i)+0.5_dp))/2.0_dp, i=1,N_Center)
!        write(FMT,'(a,i,a)')'(',N_Pair,'F10.4)'
!        write (unit_out,FMT) ( J_values(i), i=1,N_Pair )
!        close(unit_out)
       ! print the invese matrix of A=M_i*M_j
       write(unit_out,'(a)')'#  Inverse of Matrix A: (A^-1)_(k,p) [k=BS State; p=Pair]'
       write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
       call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
       write(unit_out,*)

    END SUBROUTINE Calculate_J

END MODULE J_SVD

