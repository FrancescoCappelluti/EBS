MODULE CG_TOOLS

USE PARAMETERS

IMPLICIT NONE

PRIVATE

real(kind=dp), PARAMETER, PRIVATE :: EPS = 0.0001_dp

PUBLIC :: Calculate_CG

CONTAINS





SUBROUTINE Calculate_CG(state,jj_list,CG_list,size_list)
real(kind=dp),dimension(:),intent(in)    ::   state    !|s1s2m1m2>
real(kind=dp), dimension(:),intent(out)  :: jj_list
real(kind=dp), dimension(:),intent(out)  :: CG_list
INTEGER, intent(out)            :: size_list           !dim/totalnumber of CG

real(kind=dp) :: j1,j2,m1,m2
real(kind=dp) :: mm,jj,mm_max,mm_min,jj_min,jj_max
real(kind=dp) :: a,b,c,d,z,Spin_SQR
INTEGER :: i,k,n_cg,N_jj,N_mm

!! -------------------------------------------------- !!
!             s1,s2,s3,...,sn;m1,m2,m3,....mn
!           ASSIGNING VALUES FROM THE BASIS SET
!! -------------------------------------------------- !!
  j1=state(1)
  j2=state(2)
  m1=state(3)
  m2=state(4)

  a=0.0_dp;b=0.0_dp;c=0.0_dp;d=0.0_dp;z=0.0_dp;Spin_SQR=0.0_dp
  CG_list(:)=0.0_dp
  jj_list(:)=0.0_dp
  n_cg=0

!! -------------------------------------------------- !!
!     CONDITION FOR RANGE OF SPINS AND BASIS SET
!         m = |s1 - s2| , s1 + s2
!! -------------------------------------------------- !! 

  jj_min = abs(j1-j2)
  jj_max = j1+j2
  N_jj= int((jj_max - jj_min)+ 1)
  do i=0,N_jj-1
    jj=i+jj_min
    mm_min=(-jj)
    mm_max=(+jj)
    N_mm=int(mm_max-mm_min+1)
    do k=0,N_mm-1
      mm=k+mm_min
      d=0.0_dp

!! -------------------------------------------------- !!
!   OBTAINING THE LISTS OF RESULTS FOR CLEBSCH-GORDAN
!        COEFFIECENT AND TOTAL SPIN, S
!! -------------------------------------------------- !! 

      if(CHECK(m1,m2,j1,j2,mm,jj))then
          CALL Calculate_one_CG(m1,m2,j1,j2,mm,jj,d)
          n_cg=n_cg+1
          CG_list(n_cg)=d
          jj_list(n_cg)=real(jj,kind=dp)
       end if
    end do
  end do
  size_list=n_cg


END SUBROUTINE Calculate_CG



SUBROUTINE Calculate_one_CG(v_m1,v_m2,v_j1,v_j2,v_mm,v_jj,d)
real(kind=dp),INTENT(IN)   :: v_m1,v_m2,v_j1,v_j2
real(kind=dp),INTENT(IN)  :: v_jj,v_mm
real(kind=dp)              :: m1,m2,mm,j1,j2,jj
real(kind=dp)              :: segno,tmp,a,b,c
real(kind=dp),INTENT(OUT)  :: d
INTEGER :: k_min,k_max,z
     m1=v_m1
     m2=v_m2
     mm=v_mm
     j1=v_j1
     j2=v_j2
     jj=v_jj
     segno=1.0_dp

!! -------------------------------------------------- !!
!  CONDITIONS FOR CALCULATING CLEBSCH-GORDAN COEFFIECENT
!     FOR SOLUTIONS with m < 0 and j1 < j2
!! -------------------------------------------------- !!

     if (abs(m1+m2-mm)<eps) then
       if(mm.lt.-eps)then
         m1=-v_m1
         m2=-v_m2
         mm=-v_mm
         if(abs(mod(v_jj-v_j1-v_j2,2.0_dp)).gt.eps)then
           segno=-1.0_dp
         endif
       endif
       if(j2-j1.gt.eps)then
          j1=v_j2
          j2=v_j1
          m1=v_m2
          m2=v_m1
          if(abs(mod(v_jj-v_j1-v_j2,2.0_dp)).gt.eps)then
             segno=-1.0_dp
           endif
       endif

!! -------------------------------------------------- !!
! CONDITIONS FOR CALCULATING CLEBSCH-GORDAN COEFFIECENT
!    UNDER THE SUMMATION SIGN
!! -------------------------------------------------- !!

       k_min=0
       if((-jj+j2-m1).gt.k_min)then
         k_min=int(-jj+j2-m1)
       endif
       if((-jj+j1+m2).gt.k_min)then
         k_min=int(-jj+j1+m2)
       endif

       k_max=int(j1+j2-jj)
       if((j1-m1).lt.k_max)then
          k_max=int(j1-m1)
       endif
       if((j2+m2).lt.k_max)then
          k_max=int(j2+m2)
       endif

!! -------------------------------------------------- !!
!     CLEBSCH-GORDAN COEFFIECENT FORMULATION
!! -------------------------------------------------- !!

       a=sqrt(  ( (2.0_dp*jj+1.0_dp)*R_Fact(jj+j1-j2)*R_Fact(jj-j1+j2)*R_Fact(j1+j2-jj) ) / (R_Fact(j1+j2+jj+1.0_dp) ) )
       b=sqrt(R_Fact(jj+mm)*R_Fact(jj-mm)*R_Fact(j1-m1)*R_Fact(j1+m1)*R_Fact(j2-m2)*R_Fact(j2+m2))
       c=0.0_dp
       do z=k_min,k_max
          tmp=(1.0_dp/(R_Fact(real(z,kind=dp))*R_Fact(j1+j2-jj-z)*R_Fact(j1-m1-z)*&
            R_Fact(j2+m2-z)*R_Fact(jj-j2+m1+z)*R_Fact(jj-j1-m2+z)))
          if(mod(z,2)>eps)then
            tmp=-tmp
          endif
          c=c+tmp
        enddo
        d = a*b*c

     endif
        d=segno*d

END SUBROUTINE Calculate_one_CG


!! --------------------------------------------------- !!
! CONDITIONS FOR THE NON-OCCURANCE OF CLEBSCH-GORDAN 
!    COEFFIECENT
!! --------------------------------------------------- !!

FUNCTION CHECK(m1,m2,j1,j2,mm,jj) RESULT(test)
real(kind=dp)    :: m1,m2,mm,j1,j2,jj
LOGICAL :: test

   test=.TRUE.
   if(abs(mm)>jj) then
        test=.FALSE.
   elseif (MOD((j1-j2-mm),1.0_dp) >eps) then
        test=.FALSE.
   elseif (MOD((j1-j2+mm),1.0_dp) >eps) then
        test=.FALSE.
   elseif (abs(m1+m2-mm) > eps) then
        test=.FALSE.
   endif

END FUNCTION CHECK


!! --------------------Factorial--------------------- !!
!    Function to calculate factorials resursively
!! -------------------------------------------------- !!

RECURSIVE FUNCTION R_Fact(R)  RESULT(Fact)
real(kind=dp), INTENT(IN) :: R
real(kind=dp) :: Fact


IF (R < 1.0_dp) THEN
   Fact = 1.0_dp
ELSE
   Fact = R * R_Fact(R-1.0_dp)
END IF

END FUNCTION R_Fact




END MODULE CG_TOOLS

