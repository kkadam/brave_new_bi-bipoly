subroutine compute_virial_error(psi, rho_1d, rho_2e, h, omega, volume_factor, virial_error1, virial_error2, virial_error, K_part, Pi_part, W_part)
implicit none
include 'runscf.h'
!include 'mpif.h'
!*****************************************************************************************
!
!  subroutine arguments
!

real, dimension(numr, numphi), intent(in) :: psi

real, dimension(numr,numz,numphi) :: h

real, intent(in) :: omega, rho_1d, rho_2e

real, intent(in) :: volume_factor

real, intent(out) :: virial_error1

real, intent(out) :: virial_error2

real, intent(out) :: virial_error

real, intent(out) ::  K_part, Pi_part, W_part

!
!*****************************************************************************************
!
!  global variables
!

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr) :: rhf, r, rhfinv, rinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

!
!*****************************************************************************************
!
! local variables
!

real, dimension(numr,numz,numphi) :: temp

real :: w1,  w2, wtot

real :: t1, t2, ttot

real :: s1, s2, stot

real :: ret1, ret2, global_ret1, global_ret2

integer :: I, J, K

!
!*****************************************************************************************

virial_error = 0.0
virial_error1 = 0.0
virial_error2 = 0.0

! sum up the virial pressuure
       do i = phi2, phi3
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_2e) then
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(nc2+1.0)
               else
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(n2+1.0)
               endif
             enddo
          enddo
       enddo
       do i = phi4, numphi
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_1d) then
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(nc1+1.0)
               else
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(n1+1.0)
               endif
             enddo
          enddo
       enddo
       do i = 1, phi1
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_1d) then
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(nc1+1.0)
               else
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(n1+1.0)
               endif
             enddo
          enddo
       enddo

call binary_sum(temp, ret1, ret2)
s1 = volume_factor * ret1 
s2 = volume_factor * ret2 
stot = s1 + s2

! sum up the potenntial energy
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * pot(I,J,K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
w1 = 0.5 * volume_factor * ret1
w2 = 0.5 * volume_factor * ret2
wtot = w1 + w2

! sum up the  kinetic ennergy of rotation
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * psi(I,K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
t1 = - omega * omega * volume_factor * ret1
t2 = - omega * omega * volume_factor * ret2
ttot = t1 + t2

virial_error  = abs(2.0*ttot + 3.0*stot + wtot) / abs(wtot)
virial_error1 = abs(2.0*t1   + 3.0*s1   + w1  ) / abs(w1)
virial_error2 = abs(2.0*t2   + 3.0*s2   + w2  ) / abs(w2)

K_part = t1 + t2
Pi_part = s1 + s2
W_part = w1 + w2

end subroutine compute_virial_error
