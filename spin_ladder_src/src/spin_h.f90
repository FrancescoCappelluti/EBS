MODULE SPIN_H

USE PARAMETERS, ONLY: dp
USE CG_TOOLS,   ONLY: Calculate_CG

IMPLICIT NONE

PRIVATE

real(kind=dp), PARAMETER, PRIVATE :: EPS = 0.0001_dp

PUBLIC :: get_Basis_Set,&
          H_calculate_one_element,&
          get_spin_ladder,&
          S2_calculate_one_element

CONTAINS

  SUBROUTINE H_calculate_one_element(state_l,state_r,spin_array,J_list,M_jj,M_cg,H_value)
real(kind=dp), DIMENSION(:), intent(in) :: state_l,state_r!Fe4S4 case: one of state_l=[-2.00 -2.50 -2.50 -2.00]
real(kind=dp), DIMENSION(:), intent(in) :: spin_array     !Fe4S4 case: spin_array=[2,5/2,5/2,2]
real(kind=dp), DIMENSION(:), intent(in) :: J_list
real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_jj
real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_cg
real(kind=dp),intent(out)               :: H_value

real(kind=dp),dimension(size(J_list))    :: H_pairs
real(kind=dp), dimension(4)            :: pair_l,pair_r   !pair_l=<s1s2m1m2|; pair_r = |s1s2m1m2>; |2,5/2,-2,-5/2> for Fe4S4 case
real(kind=dp)                          :: Delta_Dirac
!real(kind=dp), DIMENSION(:,:,:,:,:)    :: M_CG
INTEGER :: num_spin,num_pair
INTEGER :: i,j,k,k_pair

num_spin=size(spin_array)
num_pair=size(J_list)
H_pairs(:)=0.0_dp

!! ------ Decompose Hamiltonian in Spin Pair terms ----- !!
!  < L |- J*CG*[(Stot^2 - s1^2 -s2^2)]| R >
!  < L |- J*CG*[(Stot^2 - s1^2 -s3^2)]| R > ...
!  < L |- J*CG*[(Stot^2 - sn^2 -sn+1^2)]| R >             !
!          Sum each term to obtain the result
!! ----------------------------------------------------- !!

!! ----------------------------------------------------- !!
!        CONDITION FOR SUMMATION OF THE HALMITONIAN
! <m3_l|m3_r>*...*<mn_l|mn_r>*<m1_l m2_l|- J*CG*[(Stot^2 - s1^2 -s2^2)]|m1_r m2_r>
! <m2_l|m2_r>*...*<mn_l|mn_r>*<m1_l m3_l|- J*CG*[(Stot^2 - s1^2 -s3^2)]|m1_r m3_r> ...
! <m1_l|m1_r>*...*<mn-1_l|mn-1_r>*<m_nl mn+1_l|- J*CG*[(Stot^2 - sn^2 -sn+1^2)]|mn_r mn+1_r>
! Therefore  for the <L| or |R> if m3_l = m3_r implies <m3_l|m3_r> = 1 and 0 otherwise (Dirac delta)
!! ------------------------------------------------------  !!


k_pair=0
do i=1,num_spin
   do j=i+1,num_spin
      k_pair=k_pair+1
      !Apply Orthogonality to the states not involved in this pair
      Delta_Dirac=1.0_dp
      do k=1,num_spin
         if((k.ne.i).and.(k.ne.j))then
            if( .not.( abs( state_l(k)-state_r(k) ).lt.EPS ) )Delta_Dirac=0.0_dp
         endif
      enddo
      if(Delta_Dirac.gt.EPS)then
         !S1 and S2 for a specific pair
         pair_l(1)=spin_array(i)
         pair_r(1)=spin_array(i)
         pair_l(2)=spin_array(j)
         pair_r(2)=spin_array(j)
         !M1 and M2 for that specif pair
         pair_l(3)=state_l(i)
         pair_r(3)=state_r(i)
         pair_l(4)=state_l(j)
         pair_r(4)=state_r(j)

         Call calculate_spinpair_coupling(i,j,pair_l,pair_r,J_list(k_pair),M_jj,M_cg,H_pairs(k_pair))
      endif
   enddo
enddo

H_value=0.0_dp
do k_pair=1,num_pair
 H_value = H_value + H_pairs(k_pair)
enddo



CONTAINS

  SUBROUTINE calculate_spinpair_coupling(s_idx,s_jdx,pair_l,pair_r,J_coupling,M_jj,M_cg,val_out)
    Integer, intent(in) :: s_jdx,s_idx
    real(kind=dp), intent(in) :: J_coupling
    real(kind=dp), DIMENSION(:), intent(in) :: pair_l,pair_r  !(2,2.5,-2,-2.5)
    real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_jj
    real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_cg
    real(kind=dp), intent(out)   :: val_out
    real(kind=dp) :: jj_pair_l,CG_pair_l,jj_pair_r,CG_pair_r

    Integer :: size_pair_l,size_pair_r
    Integer :: cg_i,cg_j,m_l_idx, m_l_jdx, m_r_idx ,m_r_jdx
    real(kind=dp)    :: JJ_Square,S1_Square,S2_Square
    real(kind=dp)    :: M_l,M_r

    val_out=0.0_dp


!! ----------------------------------------------------- !!
!                 < state_l| H |state_r>                   !
!  <s1 m1,s2 m2, ... ,sn mn| H |s1 m1,s2 m2, ... ,sn mn>
! H = -2 J*[(s1s2)+(s1s3)+ ... +(s1sn)+(s2s3)+ ... +(snsn+1)] !
! H = - J*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] !
!  <s1 m1,s2 m2, ... ,sn mn| - J*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] |s1 m1,s2 m2, ... ,sn mn>  !
!! ----------------------------------------------------- !!


!! ----------------------------------------------------- !!
!  CONDITIONS FOR THE HAMILTONIAN WITH THE COMBINATION OF !
!      THE STATES (<L|R>) = <S_l,M_l|S_r,M_r>             !
! IF  M_l = M_r and S_l = S_r  then <S_l,M_l|S_r,M_r> = 1 !
!! ----------------------------------------------------- !!


M_l=pair_l(3)+pair_l(4) ! 4.5
M_r=pair_r(3)+pair_r(4) ! 4.5

m_l_idx=int(pair_l(3)+pair_l(1))+1 !index of <m_i|, index of <m| = m_i + s_i + 1
m_l_jdx=int(pair_l(4)+pair_l(2))+1 !index of <m_j|
m_r_idx=int(pair_r(3)+pair_r(1))+1 !index of |m_i>
m_r_jdx=int(pair_r(4)+pair_r(2))+1 !index of |m_j>


!M_l =M_r
if(abs(M_l-M_r).lt.EPS)then
   size_pair_r=size(M_jj(s_idx,s_jdx,m_r_idx,m_r_jdx,:))
   size_pair_l=size_pair_r
   do cg_i=1,size_pair_l
      do cg_j=1,size_pair_r
         jj_pair_l=M_jj(s_idx,s_jdx,m_l_idx,m_l_jdx,cg_i)
         jj_pair_r=M_jj(s_idx,s_jdx,m_r_idx,m_r_jdx,cg_j)
         CG_pair_l=M_cg(s_idx,s_jdx,m_l_idx,m_l_jdx,cg_i)
         CG_pair_r=M_cg(s_idx,s_jdx,m_r_idx,m_r_jdx,cg_j)
         !S_l =S_r
         if(abs(jj_pair_l-jj_pair_r).lt.EPS)then
            JJ_Square=jj_pair_l*(jj_pair_l+1.0_dp)
            S1_Square=pair_l(1)*(pair_l(1)+1.0_dp)
            S2_Square=pair_l(2)*(pair_l(2)+1.0_dp)
            val_out=val_out - J_coupling * CG_pair_l * CG_pair_r * (JJ_Square-S1_Square-S2_Square)
         endif
      enddo
   enddo
endif

END SUBROUTINE calculate_spinpair_coupling

END SUBROUTINE H_calculate_one_element

SUBROUTINE get_Basis_Set(spin_array,basis_set)
real(kind=dp),dimension(:), intent(in)      :: spin_array !Fe4S4 case: Spin=[2,5/2,5/2,2]
real(kind=dp),dimension(:,:),intent(out)    :: basis_set  !one of basis_set=[-2.00 -2.50 -2.50 -2.00]
integer   :: num_spin
integer    :: num_state
integer    :: s
real(kind=dp),dimension(size(spin_array))   :: m_array

num_state=size(basis_set(1,:))
num_spin=size(spin_array)

m_array(:)=-spin_array(:)

!Creating an entire Matrix to contain both the Spins repeatition and the basis set Matrix
do s=1,num_state
!    basis_set(1:num_spin,s)=spin_array
!    basis_set(num_spin+1:2*num_spin,s)=m_array
    basis_set(:,s)=m_array
    call next(m_array,spin_array)
end do


CONTAINS


recursive subroutine next(m_array,spin_array)
  real(kind=dp),dimension(:), intent(inout) :: m_array
  real(kind=dp),dimension(:), intent(in) :: spin_array
  integer :: num_spin
  num_spin=size(spin_array)
  m_array(num_spin)=m_array(num_spin)+1.0_dp
  if((num_spin.gt.1).and.(m_array(num_spin).gt.(spin_array(num_spin)+EPS)))then
     m_array(num_spin)=-1.0_dp*spin_array(num_spin)
     call next(m_array(1:num_spin-1),spin_array(1:num_spin-1))
  endif
end subroutine next

End SUBROUTINE get_Basis_Set


SUBROUTINE S2_calculate_one_element(state_l,state_r,spin_array,M_jj,M_cg,S2_tot)
real(kind=dp), DIMENSION(:), intent(in) :: state_l,state_r
real(kind=dp), DIMENSION(:), intent(in) :: spin_array
real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_jj
real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_cg
real(kind=dp),intent(out)               :: S2_tot


real(kind=dp)                          :: S2_tot_pair
real(kind=dp), dimension(4)            :: pair_l,pair_r
real(kind=dp)                          :: Delta_Dirac

INTEGER :: num_spin
INTEGER :: i,j,k

num_spin=size(spin_array)

!! TOT_S2 = SUM S2_i  + 2 * SUM (S_i*S_j)
!! TOT_S2 = SUM S2_i  +  SUM (S2_tot - S2_i - S2_j)
!! First Term Only for diagonal elements
!! Second Term Similar to the Hamiltonian
!!
!! ------ Decompose in Spin Pair terms ----- !!
!     < L | Stot^2 | R >
!  < L | Stot(Stot+1) | R >
!! ----------------------------------------------------- !!

!! ----------------------------------------------------- !!
!        CONDITION FOR SUMMATION OF THE SECOND TERM
! <m3_l|m3_r>*...*<mn_l|mn_r>*<m1_l m2_l| Stot^2 |m1_r m2_r>
! <m2_l|m2_r>*...*<mn_l|mn_r>*<m1_l m3_l| Stot(Stot+1) |m1_r m3_r> ...
! <m1_l|m1_r>*...*<mn-1_l|mn-1_r>*<m_nl mn+1_l| Stot(Stot+1) |mn_r mn+1_r>
! Therefore  for the <L| or |R> if m3_l = m3_r implies <m3_l|m3_r> = 1 and 0 otherwise (Dirac delta)
!! ------------------------------------------------------  !!

S2_tot=0.0_dp

do i=1,num_spin
   do j=i+1,num_spin
      !Apply Orthogonality to the states not involved in this pair
      Delta_Dirac=1.0_dp
      do k=1,num_spin
         if((k.ne.i).and.(k.ne.j))then
            if( .not.( abs( state_l(k)-state_r(k) ).lt.EPS ) )Delta_Dirac=0.0_dp
         endif
      enddo
      if(Delta_Dirac.gt.EPS)then
         !S1 and S2 for a specific pair
         pair_l(1)=spin_array(i)
         pair_r(1)=spin_array(i)
         pair_l(2)=spin_array(j)
         pair_r(2)=spin_array(j)
         !M1 and M2 for that specific pair
         pair_l(3)=state_l(i)
         pair_r(3)=state_r(i)
         pair_l(4)=state_l(j)
         pair_r(4)=state_r(j)
         S2_tot_pair=0.0_dp
         ! call calculate_spinpair_coupling(pair_l,pair_r,S2_tot_pair)
         call calculate_spinpair_coupling(i,j,pair_l,pair_r,M_jj,M_cg,S2_tot_pair)
         S2_tot = S2_tot + S2_tot_pair
      endif
   enddo
enddo

! IF it is a diagonal term
Delta_Dirac=1.0_dp
do i=1,num_spin
     if( .not.( abs( state_l(i)-state_r(i) ).lt.EPS ) )Delta_Dirac=0.0_dp
enddo
if(Delta_Dirac.gt.EPS)then
  do i=1,num_spin
    S2_tot=S2_tot+(spin_array(i)*(spin_array(i)+1))
  enddo
endif



CONTAINS

  SUBROUTINE calculate_spinpair_coupling(s_idx,s_jdx, pair_l,pair_r,M_jj,M_cg, S2_tot)
    Integer, intent(in) :: s_jdx,s_idx
    real(kind=dp), DIMENSION(:), intent(in) :: pair_l,pair_r
    real(kind=dp), intent(out)   :: S2_tot

    Integer :: size_pair_l,size_pair_r
    real(kind=dp)    :: jj_pair_l,CG_pair_l,jj_pair_r,CG_pair_r
    real(kind=dp)    :: M_l,M_r
    real(kind=dp)    :: S1_Square,S2_Square,JJ_Square
    real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_jj
    real(kind=dp), DIMENSION(:,:,:,:,:), intent(in) :: M_cg
    Integer :: cg_i,cg_j, m_l_idx, m_l_jdx, m_r_idx ,m_r_jdx


    S2_tot=0.0_dp


!! ----------------------------------------------------- !!
!                 < state_l| Stot^2 |state_r>                   !
!  <s1 m1,s2 m2, ... ,sn mn| Stot^2 |s1 m1,s2 m2, ... ,sn mn>
!              Second Term
! Stot^2 = 2*[(s1s2)+(s1s3)+ ... +(s1sn)+(s2s3)+ ... +(snsn+1)] !
! Stot^2 = CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] !
!  <s1 m1,s2 m2, ... ,sn mn|  CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] |s1 m1,s2 m2, ... ,sn mn>  !
!! ----------------------------------------------------- !!


!! ----------------------------------------------------- !!
!  CONDITIONS FOR THE SPIN TOTAL WITH THE COMBINATION OF !
!      THE STATES (<L|R>) = <S_l,M_l|S_r,M_r>             !
! IF  M_l = M_r and S_l = S_r  then <S_l,M_l|S_r,M_r> = 1 !
!! ----------------------------------------------------- !!
M_l=pair_l(3)+pair_l(4) ! 4.5
M_r=pair_r(3)+pair_r(4) ! 4.5

m_l_idx=int(pair_l(3)+pair_l(1))+1 !index of <m_i|, index of <m| = m_i + s_i + 1
m_l_jdx=int(pair_l(4)+pair_l(2))+1 !index of <m_j|
m_r_idx=int(pair_r(3)+pair_r(1))+1 !index of |m_i>
m_r_jdx=int(pair_r(4)+pair_r(2))+1 !index of |m_j>

!M_l =M_r
if(abs(M_l-M_r).lt.EPS)then
   size_pair_r=size(M_jj(s_idx,s_jdx,m_r_idx,m_r_jdx,:))
   size_pair_l=size_pair_r
   do cg_i=1,size_pair_l
      do cg_j=1,size_pair_r
         jj_pair_l=M_jj(s_idx,s_jdx,m_l_idx,m_l_jdx,cg_i)
         jj_pair_r=M_jj(s_idx,s_jdx,m_r_idx,m_r_jdx,cg_j)
         CG_pair_l=M_cg(s_idx,s_jdx,m_l_idx,m_l_jdx,cg_i)
         CG_pair_r=M_cg(s_idx,s_jdx,m_r_idx,m_r_jdx,cg_j)
         !S_l =S_r
         if(abs(jj_pair_l-jj_pair_r).lt.EPS)then
            JJ_Square=jj_pair_l*(jj_pair_l+1.0_dp)
            S1_Square=pair_l(1)*(pair_l(1)+1.0_dp)
            S2_Square=pair_l(2)*(pair_l(2)+1.0_dp)
            S2_tot= S2_tot + CG_pair_l * CG_pair_r * (JJ_Square-S1_Square-S2_Square)
         endif
      enddo
   enddo
endif



END SUBROUTINE calculate_spinpair_coupling


END SUBROUTINE S2_calculate_one_element


!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  H(n,n) = real(kind=dp) symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: states(n,n) = orthonormal eigenvectors of a           !
!        energies(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!

SUBROUTINE get_spin_ladder(H,energies,states,n)
    real(kind=dp), dimension(n,n), intent(in) :: H
    real(kind=dp), dimension(n,n), intent(out) :: states
    real(kind=dp), dimension(n), intent(out) :: energies
    integer, intent(in) :: n

    integer l,inf
    real(kind=dp)  work(n*(3+n/2))

 l=n*(3+n/2)
 states=H
 call dsyev('V','U',n,states,n,energies,work,l,inf)


END SUBROUTINE get_spin_ladder

END MODULE SPIN_H
