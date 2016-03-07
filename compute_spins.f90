!************************************************************************
!* 
!*  DIAG_SPINS
!*
!************************************************************************
subroutine compute_spins(rho,omega,com,star1com,star2com,spin1,spin2,ang_mom)
implicit none
include 'runscf.h'
!************************************************************************
!*
!
!   diag_spins calculates the z component of spin angular momentum 
!   of each star about their center of mass relative to an inertial
!   frame of reference
!
!   the formula is:
!
!   spin_1 = I_1 * Omega
!
!      I_1 = SUM [ m * r^2 ]  
!          = SUM [ rho_i * volume_factor * ((x_i - x_com_1)^2 + (y_i)^2) ]
!          
!
!         
!
!*
!************************************************************************
!*
!*   Global variables
real, dimension(numr) :: r, rhf, rinv, rhfinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numphi) :: cosine
real, dimension(numphi) :: sine
common /trig/ cosine, sine

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv


!*
!************************************************************************
!*
!*   Local Variables

integer :: J, K, I
real, dimension(numr, numphi) :: xx, yy
real :: volume_factor, com_y

!*
!************************************************************************
!*
!*  Subroutine Arguments
real, dimension(numr,numz,numphi), intent(in) :: rho
real, intent(in) :: com, star1com, star2com
real, intent(in) ::  omega
real, intent(out) :: spin1, spin2, ang_mom


!*
!************************************************************************
!  initialize the local variables
   spin1 = 0.0
   spin2 = 0.0 
   ang_mom = 0.0

!------------------------------------------------------------------------------

   volume_factor = 2.0 * dr * dz * dphi


! Get cartesian coordinates 

   do i = 1, numphi
      do k = 1, numr
         xx(k,i) = ( rhf(k) * cosine(i) )
         yy(k,i) = ( rhf(k) * sine(i) )
      enddo
   enddo



! sum up the z component of spin angular momentum

          do i = 1,phi1+1
             do j = 1,numz
                do k = 1,numr
                   spin1 = spin1 + rhf(k) * rho(k,j,i) *          &
                           ( ( xx(k,i) - star1com )**2  +         &
                            yy(k,i)**2  )
                enddo
             enddo
          enddo

          do i = phi2,phi3+1
             do j = 1,numz
                do k = 1,numr
                   spin2 = spin2 + rhf(k) * rho(k,j,i) *          &
                           ( ( xx(k,i) - star2com )**2  +         &
                            yy(k,i)**2  )
                enddo
             enddo
          enddo

          do i = phi4,numphi
             do j = 1,numz
                do k = 1,numr
                   spin1 = spin1 + rhf(k) * rho(k,j,i) *          &
                           ( ( xx(k,i) - star1com )**2  +         &
                            yy(k,i)**2  )   
                enddo
             enddo
          enddo


! Calculate total orbital angular momentum
          do i = 1,numphi
             do j = 1,numz
                do k = 1,numr
                 ang_mom = ang_mom + rhf(k) * rho(k,j,i) *       &
                         ( ( xx(k,i) - com )**2  +   &
                            yy(k,i)**2  )
                enddo
             enddo
          enddo


!------------------------------------------------------------------------------
spin1 = volume_factor * spin1 * omega
spin2 = volume_factor * spin2 * omega
ang_mom = volume_factor * ang_mom * omega



end subroutine compute_spins
