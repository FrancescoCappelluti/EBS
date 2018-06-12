MODULE FORCE

USE PARAMETERS, ONLY: dp
USE CG_TOOLS,   ONLY: Calculate_CG

IMPLICIT NONE

PRIVATE

real(kind=dp), PARAMETER, PRIVATE :: EPS = 0.0001_dp

PUBLIC :: engradDiffGroundHigh

CONTAINS

SUBROUTINE engradDiffGroundHigh(state_l,state_r,spin_array,J_coupling_gradient,H_gradient_value)
real(kind=dp), DIMENSION(:), intent(in) :: state_l,state_r
real(kind=dp), DIMENSION(:), intent(in) :: spin_array
real(kind=dp), DIMENSION(:), intent(in) :: J_coupling_gradient
real(kind=dp), intent(out)              :: H_gradient_value

real(kind=dp), dimension(size(J_coupling_gradient))    :: H_pairs_gradient
real(kind=dp), dimension(4)            :: pair_l,pair_r
real(kind=dp)                          :: Delta_Dirac
INTEGER :: num_spin,num_pair
INTEGER :: i,j,k,k_pair

num_spin=size(spin_array)
num_pair=size(J_coupling_gradient)
H_pairs_gradient(:)=0.0_dp

!! ------------------------------------------------------  !!
!! gradient correction based on BS
!! ------ Decompose Hamiltonian in Spin Pair terms ----- !!
!  < L |- ∇J*CG*[(Stot^2 - s1^2 -s2^2)]| R >
!  < L |- ∇J*CG*[(Stot^2 - s1^2 -s3^2)]| R > ...
!  < L |- ∇J*CG*[(Stot^2 - sn^2 -sn+1^2)]| R >             !
!          Sum each term to obtain the result
!! ----------------------------------------------------- !!

!! ----------------------------------------------------- !!
!        CONDITION FOR SUMMATION OF THE HALMITONIAN
! <m3_l|m3_r>*...*<mn_l|mn_r>*<m1_l m2_l|- ∇J*CG*[(Stot^2 - s1^2 -s2^2)]|m1_r m2_r>
! <m2_l|m2_r>*...*<mn_l|mn_r>*<m1_l m3_l|- ∇J*CG*[(Stot^2 - s1^2 -s3^2)]|m1_r m3_r> ...
! <m1_l|m1_r>*...*<mn-1_l|mn-1_r>*<m_nl mn+1_l|- ∇J*CG*[(Stot^2 - sn^2 -sn+1^2)]|mn_r mn+1_r>
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

    pair_l(1)=spin_array(i)
    pair_r(1)=spin_array(i)
    pair_l(2)=spin_array(j)
    pair_r(2)=spin_array(j)

    pair_l(3)=state_l(i)
    pair_r(3)=state_r(i)
    pair_l(4)=state_l(j)
    pair_r(4)=state_r(j)

    Call calculate_spinpair_coupling(pair_l,pair_r,J_coupling_gradient(k_pair),H_pairs_gradient(k_pair))
  endif
 enddo
enddo

H_gradient_value = 0.0_dp
do k_pair=1,num_pair
 H_gradient_value = H_gradient_value + H_pairs_gradient(k_pair)
enddo



CONTAINS


SUBROUTINE calculate_spinpair_coupling(pair_l,pair_r,J_coupling_gradient,val_out)
real(kind=dp), DIMENSION(:), intent(in) :: pair_l,pair_r
real(kind=dp), intent(in)    :: J_coupling_gradient
real(kind=dp), intent(out)   :: val_out

real(kind=dp), dimension(1000)  :: jj_pair_l,CG_pair_l,jj_pair_r,CG_pair_r
Integer :: size_pair_l,size_pair_r
Integer :: i,j
real(kind=dp)    :: JSpin_Square,S1_Square,S2_Square
real(kind=dp)    :: M_l,M_r

val_out=0.0_dp

!! ----------------------------------------------------- !!
!                 < state_l| ∇H |state_r>                   !
!  <s1 m1,s2 m2, ... ,sn mn| ∇H |s1 m1,s2 m2, ... ,sn mn>
!  ∇H = -2 ∇J*[(s1s2)+(s1s3)+ ... +(s1sn)+(s2s3)+ ... +(snsn+1)] !
!  ∇H = - ∇J*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] !
! <s1 m1,s2 m2, ... ,sn mn| - ∇J*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] |s1 m1,s2 m2, ... ,sn mn>  !
! <s1 m1,s2 m2, ... ,sn mn| - (∂J/∂E*∇E_bs)*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] |s1 m1,s2 m2, ... ,sn mn>  !
! <s1 m1,s2 m2, ... ,sn mn| - (A*∇E_bs)*CG*[(Stot^2 - s1^2 -s2^2)+ ... +(Stot^2 - sn^2 -sn+1^2)] |s1 m1,s2 m2, ... ,sn mn>  !
!! ----------------------------------------------------- !!

M_l=pair_l(3)+pair_l(4)
M_r=pair_r(3)+pair_r(4)
if(abs(M_l-M_r).lt.EPS)then
 Call Calculate_CG(pair_l,jj_pair_l,CG_pair_l,size_pair_l)
 Call Calculate_CG(pair_r,jj_pair_r,CG_pair_r,size_pair_r)
 do i=1,size_pair_l
   do j=1,size_pair_r
     if(abs(jj_pair_l(i)-jj_pair_r(j)).lt.EPS)then
        JSpin_Square=jj_pair_l(i)*(jj_pair_l(i)+1.0_dp)
        S1_Square=pair_l(1)*(pair_l(1)+1.0_dp)
        S2_Square=pair_l(2)*(pair_l(2)+1.0_dp)

        val_out=val_out - J_coupling_gradient * CG_pair_l(i)*CG_pair_r(j) * (JSpin_Square-S1_Square-S2_Square)
     endif
   enddo
 enddo

endif

END SUBROUTINE calculate_spinpair_coupling

END SUBROUTINE engradDiffGroundHigh

END MODULE FORCE
