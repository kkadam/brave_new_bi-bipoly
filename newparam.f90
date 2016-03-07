!********************************************************************
!*
!*  newparam
!*
!********************************************************************
subroutine mass_param(rho, omega, rho_1d, rho_2e, mass1, mass2,   &
                      mass_c1, mass_c2, xavg1, xavg2,   &
                      com, separation, spin1, spin2, ang_mom )
implicit none
include 'runscf.h'
!********************************************************************
!*
! The function compute_masses calculates mass of each star and 
! the cores in the given density field. It also calculates center 
! of mass of each of the star and the system and the separation. 

!********************************************************************
!*
!*  Subroutine Arguemtns
   real, dimension(numr,numz,numphi), intent(in) :: rho
   real, intent(in) :: omega, rho_1d, rho_2e
   real, intent(out) :: mass1, mass2, mass_c1, mass_c2, xavg1, xavg2, &
                     com, separation, spin1, spin2, ang_mom

!********************************************************************
!*
!*  Local variables
   real, dimension(numr,numz,numphi) :: temp
   real :: volume_factor, ret1, ret2
   integer :: i, j, k

!**************************************************************************************
!
!   global varaibles

real, dimension(numr) :: r, rhf, rinv, rhfinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numphi) :: cosine
real, dimension(numphi) :: sine
common /trig/ cosine, sine

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

!********************************************************************
!*
! Initialize variables
   mass1 = 0.0
   mass2 = 0.0
   mass_c1 = 0.0
   mass_c2 = 0.0
   xavg1 = 0.0 
   xavg2 = 0.0
   com = 0.0
   separation = 0.0
   spin1 = 0.0
   spin2 = 0.0
   ang_mom =0.0
   volume_factor = 2.0 * dr * dz * dphi

! Calculate the total mass for each star  and core mass for each star
   temp = 0.0
   do K = philwb, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)

   mass1 = volume_factor * ret1
   mass2 = volume_factor * ret2

   temp = 0.0
   do i = 1, numphi
      do j = 2, numz
         do k = 2, numr
            if (rho(k,j,i).gt.rho_1d) then
               temp(k,j,i) = rhf(k)*rho(k,j,i)
            endif
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)
   mass_c1 = volume_factor*ret1

   temp=0.0
   do i = 1, numphi
      do j = 2, numz
         do k = 2, numr
            if (rho(k,j,i).gt.rho_2e) then
                temp(k,j,i) = rhf(k)*rho(k,j,i)
            endif
         enddo
       enddo
   enddo
   call binary_sum(temp, ret1, ret2)
   mass_c2 = volume_factor*ret2

!!!!!!!!

   do i = 1, numphi
      do j = 2, numz
         do k = 2, numr
             rho(k,j,i) * xx() * rhf(

         enddo
       enddo
   enddo






! Calculate the center of mass for each star
   temp = 0.0
   do K = philwb, phiupb
      do J  = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rhf(I) * cosine(K) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)

   xavg1 = volume_factor * ret1 / mass1
   xavg2 = volume_factor * ret2 / mass2
   separation = xavg1 - xavg2
   com = separation * mass2 / (mass1 + mass2 )
   com = xavg1 - com

call compute_spins(rho,omega,com,xavg1,xavg2,spin1,spin2,ang_mom)


end subroutine newparam
