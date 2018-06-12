!/************************************************************************
PROGRAM MAIN

  USE PARAMETERS,  ONLY: dp,STDOUT,STDERR
  USE J_SVD,       ONLY: Calculate_J
!  USE Spin_Ladder, ONLY: Calculate_Spin_Ladder
  USE SPIN_H,      ONLY: get_Basis_Set,&
                         H_calculate_one_element,&
                         get_spin_ladder,&
                         S2_calculate_one_element
  USE PRINT_TOOLS, ONLY: write_matrix
  USE FORCE,       ONLY: engradDiffGroundHigh


        IMPLICIT NONE
!************************************************************************
        real(kind=dp), ALLOCATABLE :: energies(:)
        real(kind=dp), ALLOCATABLE :: M_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: A_inverse_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: M_ave(:),Spin(:)
        real(kind=dp), ALLOCATABLE :: J_values(:)
        real(kind=dp), ALLOCATABLE :: Basis_Set(:,:)
        real(kind=dp), ALLOCATABLE :: H_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: Eigen_States(:,:),Eigen_Values(:)
        real(kind=dp), ALLOCATABLE :: S2_matrix(:,:),S2_diag(:,:),S2_values(:)
        real(kind=dp), ALLOCATABLE :: Multiplicity(:)
        real(kind=dp), ALLOCATABLE :: State_l(:),State_r(:)
        real(kind=dp), ALLOCATABLE :: J_coupling_engrad(:)
        real(kind=dp), ALLOCATABLE :: engrad_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: GS_engrad(:)
        real(kind=dp), ALLOCATABLE :: GS_HS_engrad_diff(:)
        real(kind=dp) :: GS_HS_diff
        INTEGER :: N_BS_state,N_center,N_engrad,N_pair,Hilbert_Dim

        real(kind=dp)              :: val
        INTEGER :: i,j,k,nengrad
        INTEGER   ::    iostat
        INTEGER, PARAMETER     :: unit_in=27,unit_out=29
        CHARACTER(len=20)  :: FMT


        !!!!!!!!!! READ INPUT VALUES !!!!!!!!!!!!!!!!!!!!!!!!
        !write(*,'(a)',advance='no')'Reading Inputs...                '
        ! OPEN the M_values.dat and read first line with N&M
        iostat=0
        open(unit=unit_in,file='M_values.dat',status='old',iostat=iostat)
        if(iostat>0)STOP "Error Opening File M_values.dat"
        read(unit_in,*)N_BS_state,N_center
        allocate(M_matrix(N_BS_state,N_center)); M_matrix=0.0_dp

        ! Read all M_i values
        do i=1,N_BS_state
          read(unit_in,*)(M_matrix(i,j),j=1,N_center)
        enddo
        close(unit_in)

        ! OPEN the Energies.dat and read
        allocate(energies(N_BS_state)); energies=0.0_dp
        iostat=0
        open (unit_in,file='Energies.dat',status='old',iostat=iostat)
        do i=1,N_BS_state
          read (unit_in,*) energies(i)
        enddo
        close(unit_in)

        !READ energy gradient information from engrad.dat
        iostat=0
        open(unit=unit_in,file='engrad.dat',status='old',iostat=iostat)
        if(iostat>0)STOP "Error Opening File engrad.dat"
        read(unit_in,*)N_BS_state,N_engrad
        allocate(engrad_matrix(N_BS_state,N_engrad)); engrad_matrix=0.0_dp

        do i=1,N_BS_state
           read(unit_in,*)(engrad_matrix(i,j),j=1,N_engrad)
        enddo
        close(unit_in)


        !Calculating exchange coupling constant J
        write(STDOUT,'(a)',advance='no')'Calculating J and ∂J/∂E via SVD ...'
        open(unit=unit_out, file="svd_matrixes.out", action="write", status="replace")

        N_pair=N_center*(N_center-1)/2
        allocate(J_values(N_pair)); J_values=0.0_dp
        allocate(M_ave(N_center));  M_ave=0.0_dp
        allocate(A_inverse_matrix(N_BS_state,N_pair)); A_inverse_matrix=0.0_dp

        CALL Calculate_J(energies,M_matrix,A_inverse_matrix,M_ave,J_values,N_BS_state,N_center,N_pair,unit_out)
        close(unit_out)
        allocate(Spin(N_center))
        do i=1,N_center
          Spin(i)=real(int(M_ave(i)+0.5_dp))/2.0_dp
        enddo
!        open (unit=unit_out,file="SpinLadder.inp",action="write",status="replace")
!        write (unit_out,'(I3)')N_Center
!        write(FMT,'(a,i,a)')'(',N_Center,'F5.1)'
!        write (unit_out,FMT)(Spin(i), i=1,N_Center)
!        write(FMT,'(a,i,a)')'(',N_Pair,'F10.4)'
!        write (unit_out,FMT) ( J_values(i), i=1,N_Pair )
!        close(unit_out)
        write(STDOUT,*)'Done!'

!!!!!!!!!!!! DEFINE BASIS SET !!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Defining Basis-Set ...             '
        Hilbert_Dim=1
        do i=1,N_center
            Hilbert_Dim = Hilbert_Dim * int(2.0_dp * Spin(i) + 1.0_dp)
        end do

        allocate(Basis_Set(N_center,Hilbert_Dim))
        Basis_Set(:,:)=0.0_dp

        Call get_Basis_Set(Spin,Basis_Set)

        open(unit=unit_out, file="Basis_Set.out", action="write", status="replace")

        write(unit_out,'(a)')'#       Basis Set' 
        write(unit_out,*)
        write(FMT,"(a,I2,a)") '(I3,a,',(2*N_center),'F6.2)'
        do i=1,Hilbert_Dim
            write(unit_out,FMT) i,': ',Spin(:),Basis_Set(:,i)
        end do
        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!! CALCULATE HAMILTONIAN MATRIX !!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Calculating Spin Hamiltoninan ...  '
        allocate(State_l(N_center),State_r(N_center))
        State_l(:)=0.0_dp
        State_r(:)=0.0_dp
        Allocate(H_matrix(Hilbert_Dim,Hilbert_Dim))
        H_matrix(:,:)=0.0_dp

        do i=1,Hilbert_Dim
            do j=1,Hilbert_Dim
                State_l(:)=Basis_Set(:,i)
                State_r(:)=Basis_Set(:,j)
                call H_calculate_one_element(State_l,State_r,Spin,J_values,val)
                H_matrix(j,i)=val
            end do
        end do
        open(unit=unit_out, file="SpinHamiltonian.out", action="write", status="replace")
        write(unit_out,'(a)')'#     HAMILTONIAN MATRIX'
        write(unit_out,*)
        call write_matrix(matrix=H_matrix,cutoff=6,unit=unit_out)
        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!! CALCULATE SPIN LADDER !!!!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Calculating Spin Ladder ...        '
        !allocate(H_diag(num_states,num_states))
        allocate(Eigen_States(Hilbert_Dim,Hilbert_Dim))
        allocate(Eigen_Values(Hilbert_Dim))
        !H_diag(:,:)=0.0_dp
        Eigen_Values(:)=0.0_dp
        Eigen_States(:,:)=0.0_dp
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! diagonalize the calculated hamitonian matrix to get the eigen values(spin ladder
        ! GS energy respect to HS)

        call get_spin_ladder(H_matrix,Eigen_Values,Eigen_States,Hilbert_Dim)
        Eigen_States=TRANSPOSE(Eigen_States)

        !!! TEST RESULTS
        !H_diag=MATMUL(Basis_Set,MATMUL(H_matrix,TRANSPOSE(Basis_Set)))
        !call write_matrix(H_diag,5)

        open(unit=unit_out, file="EigenStates.out", action="write", status="replace")
        !!! PRINT EACH EIGENVECTOR as EACH COLUMN
        write(unit_out,'(a)')'#         EIGENSTATES'
        write(unit_out,*)
        call write_matrix(matrix=Eigen_States,cutoff=6,unit=unit_out)
        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! CALCULATE S2 TOTAL !!!!!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Calculating the Spin^2 ...         '
        Allocate(S2_matrix(Hilbert_Dim,Hilbert_Dim))
        !Define S2_matrix for Basis_Set
        S2_matrix(:,:)=0.0_dp
        do i=1,Hilbert_Dim
            do j=1,Hilbert_Dim
                State_L(:)=Basis_Set(:,i)
                State_R(:)=Basis_Set(:,j)
                call S2_calculate_one_element(State_L,State_R,Spin,val)
                S2_matrix(i,j)=val
            end do
        end do


        allocate(S2_diag(Hilbert_Dim,Hilbert_Dim))
        allocate(S2_values(Hilbert_Dim))
        S2_diag(:,:)=0.0_dp
        S2_values(:)=0.0_dp

        ! S2_matrix Diagonalisation rotating in the Hilbert space 
        S2_diag=MATMUL(Eigen_States,MATMUL(S2_Matrix,TRANSPOSE(Eigen_States)))
        !call write_matrix(S2_diag,5)

        do i=1,Hilbert_Dim
         S2_values(i)=S2_diag(i,i)
        enddo
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!  EIGEN ENERGY and MULTIPLICITY ASSIGNIMENT to SPINLADDER !!!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Assigning Spin ...                 '

        allocate(Multiplicity(Hilbert_Dim))
        Multiplicity(:)=0.0_dp


        Multiplicity(:) = (- 1.0_dp + sqrt(1.0_dp + 4.0_dp*S2_values(:)))/2.0_dp
        Multiplicity(:) = 2.0_dp * Multiplicity(:) + 1.0_dp
        open(unit=unit_out, file="SpinLadder.out", action="write", status="replace")

        write(unit_out,'(a)')'#STATE  ENERGY  SPIN'
        write(unit_out,*)
        deallocate(S2_diag);allocate(S2_diag(2,Hilbert_Dim))
        val=MAXVAL(ABS(Eigen_Values))
        j=0
        do while (val.gt.10)
          val=val/10.0_dp
          j=j+1
        enddo
        val=MAXVAL(ABS(Multiplicity))
        k=0
        do while (val.gt.10)
          val=val/10.0_dp
          k=k+1
        enddo
        j=j+10
        k=k+6
        write(FMT,"(a,I2,a,I2,a)") '(I3,a,F',j,'.6,F',k,'.2)'
        do i=1,Hilbert_Dim
           write(unit_out,FMT)i,': ',Eigen_Values(i),Multiplicity(i)
        enddo

        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!! CALCULATE FORCE of GROUND STATE !!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Calculating Ground State Force ... '
        allocate(J_coupling_engrad(N_pair)); J_coupling_engrad=0.0_dp
        allocate(GS_engrad(N_engrad)); GS_engrad=0.0_dp
        allocate(GS_HS_engrad_diff(N_engrad));GS_HS_engrad_diff=0.0_dp
        GS_HS_diff=0.0_dp


        ! This implementation was inefficient than the  SUM(∂J/∂E) by BS state implementation.

        do nengrad=1,N_engrad
            do k=1,N_BS_state
                  J_coupling_engrad(:)= A_inverse_matrix(k,:)
               do i=1,Hilbert_Dim
                  do j=1,Hilbert_Dim
                     State_l(:)=Basis_Set(:,i)
                     State_r(:)=Basis_Set(:,j)
                     call engradDiffGroundHigh(State_l,State_r,Spin,J_coupling_engrad(:),GS_HS_diff)
                  enddo
               enddo
               !print the contributions from different BS states
               write(STDOUT,*)'GS_engrad Contribution from Broken State',k, GS_HS_diff
                GS_HS_engrad_diff(nengrad)=GS_HS_engrad_diff(nengrad) + GS_HS_diff*engrad_matrix(k,nengrad)
            enddo
            GS_engrad(nengrad)=engrad_matrix(1,nengrad) + GS_HS_engrad_diff(nengrad)
        end do
        open(unit=unit_out, file="Force_Ground_States.out", action="write", status="replace")
        write(unit_out,'(a)')'#      High Spin State energy gradient'
!       print new engrad in 6*14 matrix
!       write(unit_out,'(6F10.8)') (engrad_matrix(1,nengrad),nengrad=1,N_engrad)
!       print new engrad in vertical line
        write(unit_out,'(F15.12)') (engrad_matrix(1,nengrad),nengrad=1,N_engrad)
        write(unit_out,'(a)')'#      Ground State energy gradient'
!       print new engrad in horizontal line
!       write(unit_out,'(84F10.8)') (GS_engrad(nengrad),nengrad=1,N_engrad)

        write(unit_out,'(F15.12)') (GS_engrad(nengrad),nengrad=1,N_engrad)
        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


write(STDOUT,'(a)')'############## Code Ended ##############'




END PROGRAM MAIN
