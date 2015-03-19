!*************************************************************************
!*
!*  OUTPUT
!*
!*************************************************************************
subroutine output(frnum, out_template, data_array)
implicit none
include 'runscf.h'
!*************************************************************************
!*
!*  Global Variables


!*
!************************************************************************* 
!*
!*   Local variables

character(len=54) :: out_file
integer :: i,j                     
!*
!*************************************************************************
!*
!*   Subroutine Arguments

integer :: frnum

character(len=50) :: out_template

real, dimension(numr,numz,numphi) :: data_array

!*
!*************************************************************************
! Initialize local variables
!out_template   = 'frame'

! create the filenames for the files every pe is going
write(out_file,'(a,i4)') trim(out_template),frnum 
open(unit=50,file=trim(out_file),form='unformatted',convert='BIG_ENDIAN',status='unknown') 
write(50) data_array
close(50)


!open(unit=50,file="den",form='unformatted',status='unknown')
!write(50) data_array
!close(50)

  open(unit=10,file="star1")
    do j=1,numz
      do i=1,numr
        write(10,*) i,j,data_array(i,j,1)
      enddo
      write(10,*)
    enddo
  close(10)
  print*, "File star1 printed"

!Vertical cross section of star2 (gnuplot)
  open(unit=10,file="star2")
    do j=1,numz
      do i=1,numr
        write(10,*) i,j,data_array(i,j,numphi/2)
      enddo
      write(10,*)
    enddo
  close(10)
  print*, "File star2 printed"

return
end subroutine output
