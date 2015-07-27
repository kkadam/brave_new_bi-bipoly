subroutine compute_pressure(rho,pres,kappac1,kappae1,kappac2,kappae2,rho_1d,rho_c1d,rho_2e,rho_c2e)
implicit none
include 'runscf.h'
  real, dimension(numr,numz,numphi) :: rho, pres
  real :: rho_c1d, rho_1d, rho_c2e, rho_2e
  real :: gammac1, gammac2, gammae1,gammae2
  real :: kappac1,kappac2, kappae1,kappae2
  real :: rhoth1, rhoth2
  integer :: i,j,k

!
!**************************************************************************************************
!

   gammae1 = 1.0 + 1.0/n1
   gammae2 = 1.0 + 1.0/n2
   gammac1 = 1.0 + 1.0/nc1
   gammac2 = 1.0 + 1.0/nc2


rhoth1=(rho_1d+rho_c1d)/2.0
rhoth2=(rho_2e+rho_c2e)/2.0


print*, phi1,phi2,phi3,phi4
print*,gammac1, gammac2, gammae1,gammae2
print*, kappac1,kappac2, kappae1,kappae2




!Calculate pressure

          do i = 1,phi1+1
             do j = 2,numz
                do k = 2,numr
                   if (rho(k,j,i).gt.rhoth1) then
                      pres(k,j,i) = kappac1*rho(k,j,i)**(gammac1)
                   else
                      pres(k,j,i) = kappae1*rho(k,j,i)**(gammae1)
                   endif
                enddo
             enddo
          enddo

          do i = phi2,phi3+1
             do j = 2,numz
                do k = 2,numr
                   if (rho(k,j,i).gt.rhoth2) then
                      pres(k,j,i) = kappac2*rho(k,j,i)**(gammac2)
                   else
                      pres(k,j,i) = kappae2*rho(k,j,i)**(gammae2)
                   endif
                enddo
             enddo
          enddo

          do i = phi4,numphi
             do j = 2,numz
                do k = 2,numr
                   if (rho(k,j,i).gt.rhoth1) then
                      pres(k,j,i) = kappac1*rho(k,j,i)**(gammac1)
                   else
                      pres(k,j,i) = kappae1*rho(k,j,i)**(gammae1)
                   endif
                enddo
             enddo
          enddo


!Set pressure floor
          do i = 1,numphi
             do j = 2,numz
                do k = 2,numr
                   if (rho(k,j,i).lt. epsilon) then
                      pres(k,j,i)=0.0
                   endif
                enddo
             enddo
          enddo

!call output('pressure.bin',pres)


end subroutine compute_pressure
