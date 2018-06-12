!/************************************************************************
PROGRAM MAIN

  USE PARAMETERS,  ONLY: dp,STDIN,STDOUT,STDERR,fatt_au2cm1
  USE J_SVD,       ONLY: Calculate_J_Ising, Calculate_J_Heisenberg
  USE SPIN_H,      ONLY: get_Basis_Set,&
                         H_calculate_one_element,&
                         get_spin_ladder,&
                         S2_calculate_one_element
  USE PRINT_TOOLS, ONLY: write_matrix
  USE FORCE,       ONLY: engradDiffGroundHigh
  USE INPUT,       ONLY: readinput
  USE CG_TOOLS,    ONLY: Calculate_CG



        IMPLICIT NONE
!************************************************************************
        real(kind=dp), ALLOCATABLE :: energies(:)
        real(kind=dp), ALLOCATABLE :: M_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: A_inverse_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: M_ave(:),Spin(:)      !Fe4S4 case: Spin=[2,5/2,5/2,2]
        real(kind=dp), ALLOCATABLE :: J_values(:)
        real(kind=dp), ALLOCATABLE :: Basis_Set(:,:)
        real(kind=dp), ALLOCATABLE :: H_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: H_diff_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: Eigen_States(:,:), Eigen_States_t(:,:),Eigen_Values(:),Eigen_States_GS(:),Eigen_States_HS(:)
        real(kind=dp), ALLOCATABLE :: Eigen_States_GS_Mul_H(:),Eigen_States_HS_Mul_H(:)
        real(kind=dp), ALLOCATABLE :: S2_matrix(:,:),S2_diag(:,:),right_hamitonian(:,:),S2_values(:),S2_tot(:)
        real(kind=dp), ALLOCATABLE :: Multiplicity(:)
        real(kind=dp), ALLOCATABLE :: State_l(:),State_r(:) !Fe4S4 case: State_l(:)=Basis_Set(:,1)=-2.00 -2.50 -2.50 -2.00
        real(kind=dp), ALLOCATABLE :: J_coupling_engrad(:)
        real(kind=dp), ALLOCATABLE :: engrad_matrix(:,:)
        real(kind=dp), ALLOCATABLE :: GS_engrad(:),HS_engrad(:),force_diff(:)
        real(kind=dp), ALLOCATABLE :: force_GS(:),force_HS(:),BS_contrib(:)
        real(kind=dp), ALLOCATABLE :: M_cg(:,:,:,:,:)
        real(kind=dp), ALLOCATABLE :: M_jj(:,:,:,:,:)
        real(kind=dp), ALLOCATABLE :: pair_r(:)
        real(kind=dp), ALLOCATABLE :: jj_pair_r(:),CG_pair_r(:)
        Integer :: size_pair_r
        INTEGER :: N_BS_state,N_center,N_engrad,N_pair,Hilbert_Dim,N_m,Index_GS,Index_HS

        real(kind=dp)              :: val,GS_En
        INTEGER :: i,j,k,p,nengrad,Spin_Center_engrad_index,mi,mj,cg,jj
        INTEGER   ::    iostat
        INTEGER, PARAMETER     :: unit_in=27,unit_out=29
        CHARACTER(len=20)  :: FMT
        CHARACTER(len=10) :: spin_ladder_type



        !!!!!!!!!! READ INPUT VALUES !!!!!!!!!!!!!!!!!!!!!!!!
        ! OPEN the M_values.dat and read first line with N&M
        iostat=0
        open(unit=unit_in,file='M_values.dat',status='old',iostat=iostat)
        if(iostat>0)STOP "Error Opening File M_values.dat"
        read(unit_in,*)N_BS_state,N_center
        allocate(M_matrix(N_BS_state,N_center)); M_matrix=0.0_dp

        ! Read all M_i values
        do k=1,N_BS_state
          read(unit_in,*)(M_matrix(k,j),j=1,N_center)
        enddo
        close(unit_in)

        ! Read all S2_tot values
        open(unit=unit_in,file='S2_tot.dat',status='old',iostat=iostat)
        if(iostat>0)STOP "Error Opening File S2_tot.dat"
        allocate(S2_tot(N_BS_state)); S2_tot=0.0_dp
        do k=1,N_BS_state
           read(unit_in,*) S2_tot(k)
        enddo
        close(unit_in)

        !reading input values from stdin
        write(*,'(a)',advance='no')'Reading Inputs...                '
        allocate(Spin(N_center)); Spin=0.0_dp
        call readinput(Spin,spin_ladder_type)

        ! OPEN the Energies.dat and read
        allocate(energies(N_BS_state)); energies=0.0_dp
        iostat=0
        open (unit_in,file='Energies.dat',status='old',iostat=iostat)
        do k=1,N_BS_state
          read (unit_in,*) energies(k)
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
        open(unit=unit_out, file="J.out", action="write", status="replace")

        N_pair=N_center*(N_center-1)/2
        allocate(J_values(N_pair)); J_values=0.0_dp
        allocate(M_ave(N_center));  M_ave=0.0_dp
        allocate(A_inverse_matrix(N_BS_state,N_pair)); A_inverse_matrix=0.0_dp
        if(spin_ladder_type .EQ. "Ising")then
           CALL Calculate_J_Ising(energies,M_matrix,A_inverse_matrix,M_ave,J_values,N_BS_state,N_center,N_pair,unit_out)
           !estimation of the spin moment by M_values
           do i=1,N_center
              Spin(i)=real(int(M_ave(i)+0.5_dp))/2.0_dp
           enddo
        else

           ! write(STDOUT,*)'Number of spin Centers:', N_center
           if(size(Spin) .NE. N_center)STOP "The number of spin moments from inputfile conflicts with that of spin centers!"
           CALL Calculate_J_Heisenberg(energies,S2_tot,M_matrix,A_inverse_matrix,M_ave,J_values,N_BS_state,N_center,N_pair,unit_out)
           do i=1,N_center
              Spin(i)=Spin(i)/2.0_dp
           enddo
        endif
        write(STDOUT,*)'Spin Moments:', Spin
        close(unit_out)
        write(STDOUT,*)'Done!'


        !Calculating CG coefficients

        N_m=int(2*maxval(Spin)+1)
        allocate(M_cg(N_center,N_center,N_m,N_m,N_m)); M_cg=0.0_dp
        allocate(M_jj(N_center,N_center,N_m,N_m,N_m)); M_jj=0.0_dp
        allocate(jj_pair_r(N_m)); jj_pair_r=0.0_dp
        allocate(cg_pair_r(N_m)); cg_pair_r=0.0_dp
        ! allocate(M_cg(N_center,N_center,N_center,N_center,2)); M_cg=0.0_dp
        allocate(pair_r(4)); pair_r=0.0_dp
        do i=1,N_center
           do j=1,N_center
              do mi=1,int(2*Spin(i))+1
                 do mj=1,int(2*Spin(j))+1
                    pair_r(1)=Spin(i)
                    pair_r(2)=Spin(j)
                    pair_r(3)=mi-1-Spin(i)
                    pair_r(4)=mj-1-Spin(j)
                    Call Calculate_CG(pair_r,jj_pair_r,CG_pair_r,size_pair_r)
                    do jj=1,size_pair_r
                       M_jj(i,j,mi,mj,jj)=jj_pair_r(jj)
                    enddo
                    do cg=1,size_pair_r
                       M_cg(i,j,mi,mj,cg)=cg_pair_r(cg)
                    enddo
                 enddo
              enddo
           enddo
        enddo



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

        write(unit_out,'(a)')'#       Basis Set including: Spin Array + Spin State'
        write(unit_out,*)
        write(FMT,"(a,I3,a)") '(I8,a,',(2*N_center),'F6.2)'

        !print the FMT format for preview
        !write(*,*) FMT

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
! first 4 of Basis sets for Fe4S4 case:
! State_l(:)=Basis_Set(:,1)=-2.00 -2.50 -2.50 -2.00
! State_r(:)=Basis_Set(:,1)=-2.00 -2.50 -2.50 -1.00
! State_l(:)=Basis_Set(:,1)=-2.00 -2.50 -2.50  0.00
! State_r(:)=Basis_Set(:,2)=-2.00 -2.50 -2.50  1.00
!                ...
                call H_calculate_one_element(State_l,State_r,Spin,J_values,M_jj,M_cg,val)  !Fe4S4 case: Spin=[2,5/2,5/2,2]
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
        allocate(Eigen_States_t(Hilbert_Dim,Hilbert_Dim))
        allocate(Eigen_Values(Hilbert_Dim))
        allocate(Eigen_States_GS(Hilbert_Dim))
        allocate(Eigen_States_HS(Hilbert_Dim))
        allocate(Eigen_States_HS_Mul_H(Hilbert_Dim))
        allocate(Eigen_States_GS_Mul_H(Hilbert_Dim))
        allocate(HS_engrad(N_engrad))
        allocate(force_GS(N_BS_state))
        allocate(force_HS(N_BS_state))
        allocate(BS_contrib(N_BS_state))
        !H_diag(:,:)=0.0_dp
        Eigen_Values(:)=0.0_dp
        Eigen_States(:,:)=0.0_dp
        Eigen_States_t(:,:)=0.0_dp

        call get_spin_ladder(H_matrix,Eigen_Values,Eigen_States,Hilbert_Dim)
        Eigen_States_t=TRANSPOSE(Eigen_States)
        write(STDOUT,*) (Eigen_Values(i),i=1,Hilbert_Dim)
        !!! TEST RESULTS
        !H_diag=MATMUL(Basis_Set,MATMUL(H_matrix,TRANSPOSE(Basis_Set)))
        !call write_matrix(H_diag,5)

        open(unit=unit_out, file="EigenStates.out", action="write", status="replace")
        !!! PRINT EACH EIGENVECTOR as EACH COLUMN
        write(unit_out,'(a)')'#         EIGENSTATES'
        write(unit_out,*)
        call write_matrix(matrix=Eigen_States_t,cutoff=6,unit=unit_out)
        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! CALCULATE S2 TOTAL !!!!!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Calculating the Spin^2 ...         '
        allocate(S2_matrix(Hilbert_Dim,Hilbert_Dim))
        !Define S2_matrix for Basis_Set
        S2_matrix(:,:)=0.0_dp
        do i=1,Hilbert_Dim
            do j=1,Hilbert_Dim
               State_L(:)=Basis_Set(:,i)
               State_R(:)=Basis_Set(:,j)
               call S2_calculate_one_element(State_L,State_R,Spin,M_jj,M_cg,val)
               S2_matrix(i,j)=val
            end do
        end do


        allocate(S2_diag(Hilbert_Dim,Hilbert_Dim));S2_diag(:,:)=0.0_dp
        allocate(S2_values(Hilbert_Dim)) ; S2_values(:)=0.0_dp
        allocate(right_hamitonian(Hilbert_Dim,Hilbert_Dim)); right_hamitonian(:,:)=0.0_dp
        ! S2_matrix Diagonalisation rotating in the Hilbert space
!        right_hamitonian=MATMUL(S2_matrix,TRANSPOSE(Eigen_States))
        right_hamitonian=MATMUL(S2_matrix,Eigen_States)
        S2_diag=MATMUL(Eigen_States_t,right_hamitonian)

        ! following line raise segment fault in MKL, so use above section.
        !S2_diag=MATMUL(Eigen_States,MATMUL(S2_matrix,TRANSPOSE(Eigen_States)))
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

        write(unit_out,'(a)')'#STATE  ENERGY(cm^{-1})  SPIN'
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
        write(FMT,"(a,I3,a,I3,a)") '(I8,a,F',j,'.6,F',k,'.2)'
        !write(*,*) fmt # print (I3, a, Fj.6, Fk.2)
        do i=1,Hilbert_Dim
           write(unit_out,FMT) i,': ', Eigen_Values(i), Multiplicity(i)
        enddo

        close(unit_out)
        write(STDOUT,*)'Done!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!! CALCULATE FORCE of GROUND STATE !!!!!!!!!!!!!!!
        write(STDOUT,'(a)',advance='no')'Calculating Ground State Force ... '
        allocate(J_coupling_engrad(N_pair)); J_coupling_engrad=0.0_dp
        val=0.0_dp
        allocate(GS_engrad(N_engrad)); GS_engrad=0.0_dp
        allocate(force_diff(N_engrad)); force_diff=0.0_dp
        Allocate(H_diff_matrix(Hilbert_Dim,Hilbert_Dim))
        ! SUM(∂J/∂E) by BS state, ∂J/∂E= A from Calculate_J
        ! Assuming that state which has the lowest energy is the GS.
        Index_GS=1
        ! Find index of GS(with minimal spin moment) and HS(maximal spin moment) states
        ! If and only if there is only one state with highest spin moment(HS)
        Index_HS=MAXLOC(Multiplicity(:),DIM=1)
        GS_En= energies(1) - (Eigen_Values(Index_HS)-Eigen_Values(Index_GS))/fatt_au2cm1
        Eigen_States_GS(:)=Eigen_States_t(Index_GS,:)
        Eigen_States_HS(:)=Eigen_States_t(Index_HS,:)

        ! <s1 m1,s2 m2, ... ,sn mn| - A*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] |s1 m1,s2 m2, ... ,sn mn>*∇E_bs  !
        ! -SUM_k SUM_p SUM_ij((C^HS_i*C^HS_j-C^GS_i*C^GS_j)<b_i|S^2_p*∂J_p/∂E_k|b_j>*∇E^BS_k)
        Spin_Center_engrad_index=1
        H_diff_matrix=0.0_dp
        do k=1,N_BS_state
           do p=1,N_pair
              J_coupling_engrad(p)=A_inverse_matrix(k,p)
              !               write(STDOUT,*)'J_coupling_engrad(i)',i, J_coupling_engrad(i)
           enddo
           do i=1,Hilbert_Dim
              do j=1,Hilbert_Dim
                 State_l(:)=Basis_Set(:,i)
                 State_r(:)=Basis_Set(:,j)
                 call engradDiffGroundHigh(State_l,State_r,Spin,J_coupling_engrad(:),M_jj,M_cg,val)
                 H_diff_matrix(j,i)=val
              end do
           end do
           ! force_GS=MATMUL(Eigen_States_GS,MATMUL(H_diff_matrix,TRANSPOSE(Eigen_States_GS)))
           ! force_HS=MATMUL(Eigen_States_HS,MATMUL(H_diff_matrix,TRANSPOSE(Eigen_States_HS)))
           ! do i = 1,Hilbert_Dim
           ! Eigen_States_GS_Mul_H(i) = DOT_PRODUCT(H_diff_matrix(i,:), Eigen_States_GS)
           ! enddo

           ! call dgemm('N', 'N', Hilbert_Dim, 1, Hilbert_Dim, 1.0_dp, Eigen_States_GS, Hilbert_Dim+1, H_diff_matrix, Hilbert_Dim+1, 0.0_dp, Eigen_States_GS_Mul_H, Hilbert_Dim+1);

           CALL DGEMV( 'N', Hilbert_Dim, Hilbert_Dim, 1.0_dp, H_diff_matrix, Hilbert_Dim, Eigen_States_GS, 1, 0.0_dp, Eigen_States_GS_Mul_H, 1)
           force_GS(k)=DOT_PRODUCT(Eigen_States_GS_Mul_H,Eigen_States_GS)
           ! force_GS(k)=DDOT(Hilbert_Dim,Eigen_States_GS,1,Eigen_States_GS_Mul_H,1)

           CALL DGEMV( 'N', Hilbert_Dim, Hilbert_Dim, 1.0_dp, H_diff_matrix, Hilbert_Dim, Eigen_States_HS, 1, 0.0_dp, Eigen_States_HS_Mul_H, 1)
           force_HS(k)=DOT_PRODUCT(Eigen_States_HS_Mul_H,Eigen_States_HS)
           ! force_HS(k)=DDOT(Hilbert_Dim,Eigen_States_HS,1,Eigen_States_HS_Mul_H,1)
           BS_contrib(k)=force_HS(k)-force_GS(k)
        enddo
        HS_engrad(:) = engrad_matrix(1,:)
        ! force = negative(engrad)
        do k = 1,N_BS_state
           force_diff(:) = force_diff(:) + BS_contrib(k)*engrad_matrix(k,:)
        enddo
        GS_engrad(:)= HS_engrad(:) + force_diff(:)
        ! write Force_Ground_States.out
        open(unit=unit_out, file="Force_Ground_States.out", action="write", status="replace")
        write(unit_out,'(a)')'#      High Spin State energy gradient'

!       print new engrad in 6*14 matrix
!       write(unit_out,'(6F10.8)') (engrad_matrix(1,nengrad),nengrad=1,N_engrad)
!       print new engrad in vertical line
        write(unit_out,'(F15.12)') (HS_engrad(nengrad),nengrad=1,N_engrad)
        write(unit_out,'(a)') '#      Ground State energy gradient'
        write(unit_out,'(a)') '#      Ground State Energy:'
        write(unit_out,*) (GS_En)
!       print new engrad in horizontal line
!       write(unit_out,'(84F10.8)') (GS_engrad(nengrad),nengrad=1,N_engrad)

        write(unit_out,'(F15.12)') (GS_engrad(nengrad),nengrad=1,N_engrad)
        write(unit_out,'(a)')'#      Contributions from different BSs:'
        write(unit_out,'(F15.12)') (BS_contrib(k),k=1,N_BS_state)
        close(unit_out)

! write  GS.extcomp.out for orca geo opt
        open(unit=unit_out, file="GS.extcomp.out", action="write", status="replace")
        write(unit_out,*) (GS_En)

        !       write GS engrad into 3*28 matrix
        write(unit_out,'(3F16.12)') (GS_engrad(nengrad),nengrad=1,N_engrad)
        close(unit_out)
        write(STDOUT,*)'Done!'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


write(STDOUT,'(a)')'############## Code Ended ##############'




END PROGRAM MAIN
