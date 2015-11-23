subroutine newpressure(rho,pres,h,rho_1d,rho_2e)
implicit none
include 'runscf.h'
   real, dimension(numr,numz,numphi) :: rho, pres, h
   real ::  rho_1d, rho_2e
   real :: rhoth1, rhoth2
   integer :: i,j,k

!
!**************************************************************************************************
!

   rhoth1=rho_1d
   rhoth2=rho_2e


!Calculate pressure from density and enthalpy fields
   do i = 1,phi1+1
      do j = 1,numz
         do k = 1,numr
            if (rho(k,j,i).gt.rhoth1) then
               pres(k,j,i) = h(k,j,i)*rho(k,j,i)/(1.0+nc1)
            else
               pres(k,j,i) = h(k,j,i)*rho(k,j,i)/(1.0+n1)                   
            endif
         enddo
      enddo
   enddo

   do i = phi2,phi3+1
      do j = 1,numz
         do k = 1,numr
            if (rho(k,j,i).gt.rhoth2) then
               pres(k,j,i) = h(k,j,i)*rho(k,j,i)/(1.0+nc2)
            else
               pres(k,j,i) = h(k,j,i)*rho(k,j,i)/(1.0+n2)
            endif
         enddo
      enddo
   enddo

   do i = phi4,numphi
      do j = 1,numz
         do k = 1,numr
            if (rho(k,j,i).gt.rhoth1) then
               pres(k,j,i) = h(k,j,i)*rho(k,j,i)/(1.0+nc1)
            else
               pres(k,j,i) = h(k,j,i)*rho(k,j,i)/(1.0+n1)
            endif
         enddo
      enddo
   enddo


!Set pressure floor
   do i = 1,numphi
      do j = 1,numz
         do k = 1,numr
            if (rho(k,j,i).lt. epsilon) then
                pres(k,j,i)=0.0
            endif
         enddo
      enddo
   enddo


end subroutine newpressure
