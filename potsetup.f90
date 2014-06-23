!********************************************************************
!* 
!*  POTSETUP
!*
!********************************************************************
subroutine potsetup(have_green_funcs)
implicit none
include 'runscf.h'
include 'pot.h'
!********************************************************************
!*
!  potsetup does intitialization for the poisson solver
!  package.  In order of occurance this is what potsetup
!  does:
!
!     -> initialize trigonometric functions used by bessel
!
!     -> call tm and sm to fill in the Green function tmr and smz
!        respectively
!
!     -> setup the geometric arrays for the ADI solver
!
!*
!********************************************************************
!*
!*  Subroutine Arguments

logical, intent(inout) :: have_green_funcs

!*
!********************************************************************
!*
!*  Global Variables

real, dimension(numr,numz,numr,mmax) :: tmr
real, dimension(numr,numz,numz,mmax) :: smz
common /green_functions/ tmr, smz

real, dimension(numphi,mmax) :: bes_cos, bes_sin
common /bessel_trig/ bes_cos, bes_sin

real, dimension(numr) :: ar, cr, alphar
real, dimension(numr,numphi) :: brb
common /ADI_R_sweep/ ar, cr, alphar, brb

real :: az, cz
real, dimension(numr) :: alphaz, betaz
real, dimension(numz) :: bzb
real, dimension(numr,numphi) :: elambdazb
common /ADI_Z_sweep/ az, cz, alphaz, betaz, bzb, elambdazb

real :: gamma, piinv, four_pi
common /pot_constants/ gamma, piinv, four_pi

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor, &
                   viscosity

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g, r_g, rhfinv_g, rinv_g, zhf_g

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

!*
!********************************************************************
!*
!*  Local Variables
  
real :: mindex, lindex, drinv2, dphiinv2

real, dimension(numphi) :: elm, m1mode

real, dimension(numr) :: betar

integer :: J, K, L, M, lstop, index, mode

!*
!********************************************************************
!  initialize the local variables
mindex = 0.0
lindex = 0.0
drinv2 = 0.0
dphiinv2 = 0.0
do L = 1, numphi
   elm(L)    = 0.0
   m1mode(L) = 0.0
enddo
do J = 1, numr
   betar(J) = 0.0
enddo
lstop = 0
index = 0
mode = 0 

piinv = 1.0 / pi
four_pi = 4.0 * pi

!  setup the trig arrays used by bessel, use lindex and mindex to
!  avoid casts which are relatively expensive
mindex = 1.0
lindex = 1.0
if( isym == 3 ) then
   do M = 1, mmax
      do L = philwb, phiupb
         bes_cos(L,M) = cos(dphi*(mindex-1.0)*(2.0*lindex-1.0))
         bes_sin(L,M) = sin(dphi*(mindex-1.0)*(2.0*lindex-1.0))
         lindex = lindex + 1.0
      enddo
      lindex = 1.0
      mindex = mindex + 1.0
   enddo
else
   do M = 1, mmax
      do L = philwb, phiupb
         bes_cos(L,M) = cos(0.5*dphi*(mindex-1.0)*(2.0*lindex-1.0))
         bes_sin(L,M) = sin(0.5*dphi*(mindex-1.0)*(2.0*lindex-1.0))
         lindex = lindex + 1.0
      enddo
      lindex = 1.0
      mindex = mindex + 1.0
   enddo
endif

!  call tm and sm to fill in the Green functions for the top (and bottom)
!  potential - tmr and the Green function for the side potential - smz
if ( have_green_funcs ) then
   open(unit=20,file='tmr_array',form='unformatted',status='old')
   open(unit=21,file='smz_array',form='unformatted',status='old')
   read(20) tmr
   read(21) smz
   close(20)
   close(21)
else
   call tm(tmr)
   call sm(smz)
   open(unit=20,file='tmr_array',form='unformatted',status='new')
   open(unit=21,file='smz_array',form='unformatted',status='new')
   write(20) tmr
   write(21) smz
   close(20)
   close(21)
   have_green_funcs = .true.
endif

!  the rest of the code in potsetup sets up arrays that hold geometric
!  information for the discretization of Poisson's equation.  The
!  arrays are all used by the ADI solver
gamma = dzinv * dzinv
drinv2 = drinv * drinv 
dphiinv2 = dphiinv * dphiinv

lstop = numphi_by_two + 1
do L = 1, numphi
   if( isym == 3 ) then
      if( L <= lstop ) then
         mode = (L-1)*(L-1)
      else
         mode = (L-lstop)*(L-lstop)     
      endif
   else
      if( L <= lstop ) then
         mode = L - 1
      else
         mode = L - lstop
      endif
   endif
   m1mode(L) = (-1.0)**mode
   elm(L) = cos(mode*dphi)
enddo
      
do J = 2, numr-1

   alphar(J) = -r_g(J+1)*rhfinv_g(J)*drinv2

   betar(J) = -r_g(J)*rhfinv_g(J)*drinv2

enddo

! have to used indexed global radius array to initialize 
! alphaz and betaz because they will be used when the radial
! data is block distributed across numz_procs
index = 2
do J = 2, numr-1

   alphaz(J) = -r_g(index+1)*rhfinv_g(index)*drinv2

   betaz(J) = -r_g(index)*rhfinv_g(index)*drinv2

   index = index + 1

enddo

ar = betar

cr = alphar

az = - gamma
 
cz = - gamma

do L = 1, numphi
   do J = 3, numr-1
      brb(J,L) = 2.0*drinv2 - 2.0*(elm(L)-1.0)*dphiinv2*           &
                 rhfinv_g(J)*rhfinv_g(J)
   enddo
enddo
if( isym == 3 ) then
   do L = 1, numphi
      brb(2,L) = -alphar(2) - 2.0*betar(2) - 2.0*                  &
                  (elm(L)-1.0)*dphiinv2*rhfinv_g(2)*               &
                  rhfinv_g(2)
   enddo
else
   do L = 1, numphi
      brb(2,L) = -alphar(2) + (m1mode(L)-1.0)*betar(2) -           &
                  2.0*(elm(L)-1.0)*dphiinv2*rhfinv_g(2)*           &
                  rhfinv_g(2)
   enddo
endif

do K = 2, numz-1
   bzb(K) = 2.0*gamma
enddo
if( isym /= 1 ) bzb(2) = gamma

do L = 1, numphi
   index = 2
   do J = 2, numr-1
      elambdazb(J,L) = -2.0*drinv2 + 2.0*(elm(L)-1.0)*             &
                        dphiinv2*rhfinv_g(index)*                  &
                        rhfinv_g(index)
      index = index + 1
   enddo
enddo
if( isym == 3 ) then
   do L = 1, numphi
      elambdazb(2,L) = alphaz(2) + 2.0*betaz(2) + 2.0*             &
                       (elm(L)-1.0)*dphiinv2*rhfinv_g(2)*          &
                       rhfinv_g(2)
   enddo
endif
if( isym /= 3 ) then
   do L = 1, numphi
      elambdazb(2,L) = alphaz(2) - (m1mode(L)-1.0)*betaz(2) +      &
                       2.0*(elm(L)-1.0)*dphiinv2*rhfinv_g(2)*      &
                       rhfinv_g(2)
   enddo
endif
    
return
end subroutine potsetup
