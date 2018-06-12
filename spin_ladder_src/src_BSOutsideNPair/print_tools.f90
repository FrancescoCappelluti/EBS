MODULE print_tools

USE PARAMETERS, ONLY: dp,stdout

IMPLICIT NONE

PRIVATE

PUBLIC :: write_matrix

CONTAINS

subroutine write_matrix(matrix,cutoff,unit)
  real(kind=dp), dimension(:,:), intent(in) :: matrix
  integer, intent(in) :: cutoff
  integer, intent(in), optional :: unit

  real(kind=dp) :: max_val 
  integer :: nrows,ncolumns 
  integer :: nblocks
  character(len=100) :: FMT1,FMT2
  integer :: i,k,first,last
  integer :: unit_out

  if(PRESENT(unit))then
    unit_out=unit
  else
    unit_out=stdout
  endif

  nrows=size(matrix(1,:))
  ncolumns=size(matrix(:,1))
  nblocks=ncolumns/cutoff

  !Max Value Magnitude Determination
  max_val=MAXVAL(ABS(matrix))
  i=0
  do while (max_val.gt.10)
    max_val=max_val/10.0_dp
    i=i+1
  enddo

  i=i+10
  write(FMT1,"(a,I2,a,I2,a)") '(a,',cutoff,'I',i,')'
  write(FMT2,"(a,I2,a,I2,a)") '(I3,a,',cutoff,'F',i,'.6)'

  do k=1,nblocks
   first=1+((k-1)*cutoff)
   last=k*cutoff
   write(unit_out,FMT1)'     ',(i,i=first,last)
   do i=1,nrows
     write(unit_out,FMT2) i,': ',matrix(first:last,i)
   enddo
  enddo
  
  ncolumns=mod(ncolumns,cutoff)
  if(ncolumns.gt.0)then
    first=1+nblocks*cutoff
    last=ncolumns+nblocks*cutoff
    write(unit_out,FMT1)'     ',(i,i=first,last)
    do i=1,nrows
      write(unit_out,FMT2) i,': ',matrix(first:last,i)
    enddo
  endif

end subroutine write_matrix


END MODULE print_tools

