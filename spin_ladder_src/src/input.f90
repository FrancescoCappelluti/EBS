MODULE INPUT
  USE PARAMETERS, ONLY: dp,STDIN,STDOUT
  IMPLICIT none

CONTAINS

  SUBROUTINE readinput(Spin,spin_ladder_type)
    character(len=100) :: buffer, label
    integer :: pos
    integer :: ios = 0
    integer :: line = 0

    REAL(kind=dp), INTENT(OUT), DIMENSION(:)   :: Spin
    CHARACTER(len=10), INTENT(OUT)  :: spin_ladder_type


    do while (ios == 0)
       read(STDIN, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1
          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, ' 	')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('spin_ladder_type')
             read(buffer, *, iostat=ios) spin_ladder_type
             write(STDOUT,*)'Type of Spin Ladder:', spin_ladder_type
          case ('spin_moment')
             read(buffer, *, iostat=ios) Spin
             !write(STDOUT,*)'Spin Moments from input files:', Spin
          case default
             print *, 'Skipping invalid flag at line', line
          end select
       end if
    end do
  END SUBROUTINE
END MODULE INPUT
