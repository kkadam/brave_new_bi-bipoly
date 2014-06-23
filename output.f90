!*************************************************************************
!*
!*  OUTPUT
!*
!*************************************************************************
subroutine output(out_template, data_array)
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
                     
!*
!*************************************************************************
!*
!*   Subroutine Arguments

integer :: i,j

character(len=50) :: out_template

real, dimension(numr,numz,numphi) :: data_array

!*
!*************************************************************************

!Binary density fle    
  open(unit=10,file=trim(out_template),form='unformatted',convert='BIG_ENDIAN',&
       status='unknown') 
  write(10) data_array
  close(10)
  print*, "File ", trim(out_template), " printed"
  
  
!Vertical cross section of star1 (gnuplot)  
  open(unit=10,file="star1")
    do j=1,numz
      do i=1,numr  
        write(10,*) i,j,data_array(i,j,1) 
      enddo
      write(10,*)
    enddo
  close(10)    
  print*, "File star2 printed"       

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
