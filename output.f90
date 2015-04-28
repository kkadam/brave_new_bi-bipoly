!*************************************************************************
!*
!*  OUTPUT
!*
!*************************************************************************
subroutine output(out_template, file2d, data_array)
implicit none
include 'runscf.h'
!*************************************************************************
!*
!*  Global Variables


!*
!************************************************************************* 
!*
!*   Local variables

character(len=15) :: out_file
character(len=15) :: file1
character(len=15) :: file2
integer :: i,j                     
!*
!*************************************************************************
!*
!*   Subroutine Arguments

character(len=*), intent(in) :: out_template
character(len=*), intent(in) :: file2d

real, dimension(numr,numz,numphi), intent(in) :: data_array

!*
!*************************************************************************
! Initialize local variables
!out_template   = 'frame'

! create the filenames for the files every pe is going
write(out_file,'(a,i4)') trim(out_template)
open(unit=50,file=trim(out_file),form='unformatted',convert='BIG_ENDIAN',status='unknown') 
write(50) data_array
close(50)

print*, "File ", trim(out_file), " printed"

file1=trim(file2d)//'1'
file2=trim(file2d)//'2'


  open(unit=10,file=trim(file1))
    do j=1,numz
      do i=1,numr
        write(10,*) i,j,data_array(i,j,1)
      enddo
      write(10,*)
    enddo
  close(10)
  print*, "File ", trim(file1), " printed"

!Vertical cross section of star2 (gnuplot)
  open(unit=10,file=trim(file2))
    do j=1,numz
      do i=1,numr
        write(10,*) i,j,data_array(i,j,numphi/2+1)
      enddo
      write(10,*)
    enddo
  close(10)
  print*, "File ", trim(file2), " printed"

return
end subroutine output
