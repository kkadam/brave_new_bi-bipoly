subroutine newpressure(rho,pres, h, h_e1d, h_e2e, core_template)
implicit none
include 'runscf.h'
   real, dimension(numr,numz,numphi) :: rho, pres, h
   integer, dimension(numr, numz, numphi), intent(in) :: core_template

   real ::  h_e1d, h_e2e
   integer :: i,j,k

!
!**************************************************************************************************
!

!Calculate pressure from density and enthalpy fields
   do i = 1,phi1+1
      do j = 1,numz
         do k = 1,numr
            if ( core_template(k,j,i).eq.1 ) then
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
            if ( core_template(k,j,i).eq.1 ) then
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
            if ( core_template(k,j,i).eq.1 ) then
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
