!/************************************************************************
MODULE J_SVD
  USE PARAMETERS,   ONLY: dp,STDOUT,STDERR,fatt_au2cm1
  USE PRINT_TOOLS,  ONLY: write_matrix


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Calculate_J_Ising
  PUBLIC :: Calculate_J_Heisenberg
!************************************************************************
  CONTAINS

    SUBROUTINE Calculate_J_Ising(energies,M_matrix,A,M_ave,J_values,N_state,N_Center,N_Pair,unit_out)
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
!        INTEGER   ::    iostat
!        INTEGER, PARAMETER     :: unit_in=27,unit_out=29
!        CHARACTER(len=20)  :: FMT

        write(unit_out,'(a)')'###############################################################'
        write(unit_out,'(a)')'#          J CALCULATION WITH ISING MODEL                     #'
        write(unit_out,'(a)')'###############################################################'
        write(unit_out,*)

        write(unit_out,'(a)')'#     M MATRIX: M_(k,i) [k=BS State; i=Spin Center]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Center
        call write_matrix(matrix=TRANSPOSE(M_matrix),cutoff=N_Center,unit=unit_out)
        write(unit_out,*)

        ! Make averages along columns (absolute value)
        M_ave=0.0_dp
        do k=1,N_state
          do j=1,N_Center
            M_ave(j)=M_ave(j)+abs(M_matrix(k,j))
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
        ! lwork=6*N_Pair
        ! add minimum lwork value 40
        ! if (N_state > 4) then
                ! lwork=(N_state)**2 ; ifail=0
        ! else
                ! lwork=40
        ! endif
        !lwork >= max(3*min(m, n)+max(m, n), 5*min(m,n)) 
        ! https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/402436

        lwork=max(3*min(N_state, N_Pair)+max(N_state, N_Pair), 5*min(N_state,N_Pair))
        allocate(work(lwork))
        write(unit_out,'(a)')'#     A MATRIX: (M_i*M_j)_(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        A=A*2
!        E=-(2N_state)*J

call DGESVD('S','S', N_state, N_Pair, A, N_state, s, U, N_state, Vt, N_Pair, work, lwork, info)


        write(unit_out,'(a)')'### SVD ###'
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
        write(unit_out,'(6F10.5)') (s(p),p=1,N_Pair)
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
        do p=1,N_Pair
          Sigma(p,p)=1/s(p)
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
        do k=1,N_state
           do j=1,N_Pair
            A(k,j)=A(k,j)+U(k,p)*temp_ave(p,j)
           enddo
        enddo
        enddo


        J_values=0.0_dp
        do p=1,N_Pair
          do k=1,N_state
             !J = -A^{-1}*E^BS
             J_values(p)=J_values(p)-A(k,p)*energies(k)
          enddo
        enddo

        write(unit_out,'(a)')'#  Inverse of Matrix A for Calculating J: (A^-1)_(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        !J in cm^{-1} unit
        J_values=J_values*fatt_au2cm1

        !print J coupling values into J.out
        write(unit_out,'(a)')'#    J coupling (cm^{-1}):'
        write(unit_out,'(a,I5)')'# ',N_Pair
        write(unit_out,'(ES14.7)') (J_values(p),p=1,N_Pair)
        write(unit_out,*)

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

    END SUBROUTINE Calculate_J_Ising

    SUBROUTINE Calculate_J_Heisenberg(energies,S2_tot,M_matrix,A,M_ave,J_values,N_state,N_Center,N_Pair,unit_out)
        REAL(kind=dp), INTENT(IN), DIMENSION(:)   :: energies
        REAL(kind=dp), INTENT(IN), DIMENSION(:)   :: S2_tot
        REAL(kind=dp), ALLOCATABLE, DIMENSION(:)  :: spin_correlation
        REAL(kind=dp), INTENT(IN), DIMENSION(:,:) :: M_matrix
        REAL(kind=dp), INTENT(OUT), DIMENSION(:)  :: M_ave
        REAL(kind=dp), INTENT(OUT), DIMENSION(:)  :: J_values
        INTEGER, INTENT(IN)                       :: N_state,N_Center,N_Pair
        INTEGER, INTENT(IN)                       :: unit_out
        REAL(kind=dp), INTENT(OUT), DIMENSION(N_state,N_Pair) :: A
        ! REAL(kind=dp), DIMENSION(:)  :: spin_correlation
        REAL(kind=dp), DIMENSION(N_Pair) :: S
        REAL(kind=dp), DIMENSION(N_state,N_Pair) :: U
        REAL(kind=dp), DIMENSION(N_Pair,N_Pair) :: Vt
        REAL(kind=dp), DIMENSION(N_Pair,N_Pair) :: Sigma,temp_ave

        real(kind=dp), ALLOCATABLE ::  work(:)
!        INTEGER, ALLOCATABLE :: desca(:),descu(:),descvt(:)
        INTEGER :: ifail, i,j,k,p, info, lwork
!        INTEGER   ::    iostat
!        INTEGER, PARAMETER     :: unit_in=27,unit_out=29
!        CHARACTER(len=20)  :: FMT

        write(unit_out,'(a)')'###############################################################'
        write(unit_out,'(a)')'#          J CALCULATION WITH HEISENBERG MODEL                #'
        write(unit_out,'(a)')'###############################################################'
        write(unit_out,*)
        write(unit_out,'(a)')'#     M MATRIX: M_(k,i) [k=BS State; i=Spin Center]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Center
        call write_matrix(matrix=TRANSPOSE(M_matrix),cutoff=N_Center,unit=unit_out)
        write(unit_out,*)

        !calculate heisengberg spin pairs S_tot^2=SUM(S(i)^2)+2*S1S2+2*S1S3+2*S1S4....
        allocate(spin_correlation(N_Pair));spin_correlation=0.0_dp

        call heisenberg_spin_correlation(S2_tot,spin_correlation,M_matrix,N_state,N_Center,N_Pair,unit_out)

        ! Make averages along columns (absolute value)
        M_ave=0.0_dp
        do k=1,N_state
          do j=1,N_Center
            M_ave(j)=M_ave(j)+abs(M_matrix(k,j))
          enddo
        enddo
        M_ave(:)=M_ave(:)/N_state

        A=0.0_dp
        do k=1,N_state
          p=1
          do i=1,N_Center
            do j=i+1,N_Center
              A(k,p)=sign(1.0_dp,M_matrix(k,i))*sign(1.0_dp,M_matrix(k,j))*spin_correlation(p)
              p=p+1
            enddo
          enddo
        enddo


        U=0.0_dp
        Vt=0.0_dp
        s=0.0_dp
        ! lwork=6*N_Pair
! add minimum lwork value 40
        ! if (N_state > 4) then
                ! lwork=(N_state)**2 ; ifail=0
        ! else
                ! lwork=40
        ! endif

        !lwork >= max(3*min(m, n)+max(m, n), 5*min(m,n)) 
        ! https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/402436

        lwork=max(3*min(N_state, N_Pair)+max(N_state, N_Pair), 5*min(N_state,N_Pair))

        allocate(work(lwork))
        write(unit_out,'(a)')'#     A MATRIX for calculating J: A(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)
        A=A*2

call DGESVD('S','S', N_state, N_Pair, A, N_state, s, U, N_state, Vt, N_Pair, work, lwork, info)


        write(unit_out,'(a)')'### SVD ###'
        write(unit_out,'(a)')'#     U MATRIX for calculating J: U_(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(U),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        write(unit_out,'(a)')'#    Vt MATRIX for calculating J: Vt_(p1,p1) [p1,p2=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_Pair,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(Vt),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        write(unit_out,'(a)')'#    Sigma Diagonal Values for calculating J:'
        write(unit_out,'(a,I5)')'# ',N_Pair
        write(unit_out,'(6F10.5)') (s(p),p=1,N_Pair)
        write(unit_out,*)





!!!! J CALCULATION
!        write (stdout,*) 'Starting J_ij calculations.'


       !Inverto la matrice diagonale
       !Matrix A=M*M, then calculate the inverse matrix A in following section
       !Inverse diagonalize terms
        Sigma=0.0_dp
        do p=1,N_Pair
          Sigma(p,p)=1/s(p)
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
        do k=1,N_state
           do j=1,N_Pair
            A(k,j)=A(k,j)+U(k,p)*temp_ave(p,j)
           enddo
        enddo
       enddo

       write(unit_out,'(a)')'#  Inverse of Matrix A for calculating J: (A^-1)_(k,p) [k=BS State; p=Pair]'
       write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
       call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
       write(unit_out,*)

        J_values=0.0_dp
        do p=1,N_Pair
          do k=1,N_state
              !j = -A^{-1}*E^BS
              J_values(p)=J_values(p)-A(k,p)*energies(k)
          enddo
        enddo

        J_values=J_values*fatt_au2cm1

        !print J coupling values into J.out
        write(unit_out,'(a)')'#    J coupling (cm^{-1}):'
        write(unit_out,'(a,I5)')'# ',N_Pair
        write(unit_out,'(ES14.7)') (J_values(p),p=1,N_Pair)
        write(unit_out,*)

    END SUBROUTINE Calculate_J_Heisenberg

    SUBROUTINE heisenberg_spin_correlation(S2_tot,spin_correlation,M_matrix,N_state,N_Center,N_Pair,unit_out)

      !S2_tot = (SUM_{i=1}^{N_p} s_i)^2; total spin(including all N spin centers) square from DFT such as 'Expectation value of <S**2>' in ORCA or QMC...
      REAL(kind=dp), INTENT(IN), DIMENSION(:)   :: S2_tot
      REAL(kind=dp), INTENT(IN), DIMENSION(:,:) :: M_matrix
      !spin_correlation is two spin sites interacted with each other
      ! REAL(kind=dp), INTENT(OUT), ALLOCATABLE, DIMENSION(:)  :: spin_correlation
      REAL(kind=dp), INTENT(OUT), DIMENSION(:)  :: spin_correlation
      INTEGER, INTENT(IN)                       :: N_state,N_Center,N_Pair
      INTEGER, INTENT(IN)                       :: unit_out
      REAL(kind=dp), DIMENSION(N_state,N_Pair)  :: A
      REAL(kind=dp), DIMENSION(N_Pair)         :: S
      REAL(kind=dp), DIMENSION(N_state,N_Pair) :: U
      REAL(kind=dp), DIMENSION(N_Pair,N_Pair)  :: Vt
      REAL(kind=dp), DIMENSION(N_Pair,N_Pair)  :: Sigma,temp_ave

      real(kind=dp), ALLOCATABLE ::  work(:)
!      INTEGER, ALLOCATABLE :: desca(:),descu(:),descvt(:)
      INTEGER :: ifail, i,j,k,p, info, lwork
!      INTEGER   ::    iostat
!      INTEGER, PARAMETER     :: unit_in=27,unit_out=29
!      CHARACTER(len=20)  :: FMT

      A=0.0_dp
      do k=1,N_state
         p=1
         do i=1,N_Center
           do j=i+1,N_Center
             A(k,p)=sign(1.0_dp,M_matrix(k,i))*sign(1.0_dp,M_matrix(k,j))
             p=p+1
           enddo
         enddo
      enddo


      U=0.0_dp
      Vt=0.0_dp
      s=0.0_dp
      ! lwork=6*N_Pair
      ! add minimum lwork value 40
        ! if (N_state > 4) then
                ! lwork=(N_state)**2 ; ifail=0
        ! else
                ! lwork=40
        ! endif

      !lwork >= max(3*min(m, n)+max(m, n), 5*min(m,n)) 
      !https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/402436
      !https://software.intel.com/en-us/node/469236


      lwork=max(3*min(N_state, N_Pair)+max(N_state, N_Pair), 5*min(N_state,N_Pair))

        allocate(work(lwork))
        write(unit_out,'(a)')'#     A MATRIX for Calculating Spin Correlation: A(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        A=A*2

        call DGESVD('S','S', N_state, N_Pair, A, N_state, s, U, N_state, Vt, N_Pair, work, lwork, info)

        write(unit_out,'(a)')'### SVD ###'
        write(unit_out,'(a)')'#     U MATRIX for Calculating Spin Correlation: U_(k,p) [k=BS State; p=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(U),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        write(unit_out,'(a)')'#    Vt MATRIX for Calculating Spin Correlation: Vt_(p1,p1) [p1,p2=Pair]'
        write(unit_out,'(a,I5,a,I5)')'# ',N_Pair,' ',N_Pair
        call write_matrix(matrix=TRANSPOSE(Vt),cutoff=N_Pair,unit=unit_out)
        write(unit_out,*)

        write(unit_out,'(a)')'#    Sigma Diagonal Values for Calculating Spin Correlation:'
        write(unit_out,'(a,I5)')'# ',N_Pair
        write(unit_out,'(6F10.5)') (s(p),p=1,N_Pair)
        write(unit_out,*)


        Sigma=0.0_dp
        do p=1,N_Pair
          Sigma(p,p)=1/s(p)
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
        do k=1,N_state
           do j=1,N_Pair
            A(k,j)=A(k,j)+U(k,p)*temp_ave(p,j)
           enddo
        enddo
       enddo


       write(unit_out,'(a)')'#     Inverse A MATRIX for Calculating Spin Correlation: A(k,p) [k=BS State; p=Pair]'
       write(unit_out,'(a,I5,a,I5)')'# ',N_state,' ',N_Pair
       call write_matrix(matrix=TRANSPOSE(A),cutoff=N_Pair,unit=unit_out)
       write(unit_out,*)



       spin_correlation=0.0_dp
       do p=1,N_Pair
         do k=1,N_state
            !spin_pair=A^{-1}*S2_tot
            spin_correlation(p)=spin_correlation(p)+A(k,p)*S2_tot(k)
         enddo
       enddo

       write(unit_out,'(a)')'#    Spin Square of all BS states S2_tot:'
       write(unit_out,'(a,I5)')'# ', N_state
       write(unit_out,'(ES14.7)') (s2_tot(k),k=1,N_state)
       write(unit_out,*)

       write(unit_out,'(a)')'#    Spin Correlation:'
       write(unit_out,'(a,I5)')'# ', N_Pair
       write(unit_out,'(ES14.7)') (spin_correlation(p),p=1,N_Pair)
       write(unit_out,*)
    END SUBROUTINE heisenberg_spin_correlation

END MODULE J_SVD
