!********************************************************************************
!*
!*  setup
!*
!********************************************************************************
subroutine setup( have_green_funcs )
implicit none
include 'runscf.h'
!include 'mpif.h'
!********************************************************************************
!*
!*   Subroutine Arguments

logical, intent(in) :: have_green_funcs

!*
!********************************************************************************
!*
!*   Global variables

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numr) :: rhf, r, rhfinv, rinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numphi) :: cosine, sine
common /trig/ cosine, sine

integer, dimension(3) :: boundary_condition
common /boundary_conditions/ boundary_condition

!*
!*********************************************************************************
!*
!*  Local variables

integer :: J, K, L

real :: x, xinv, offset

!*
!*********************************************************************************      
!  initialize the local variables
x = 0.0
xinv = 0.0
offset = 0.0

if( isym == 2 .or. isym == 3 ) then
   boundary_condition(1) = 1
endif

!  setup coordinates, all pe's use their own portion of the domain
!  decomposition except for the global ( _g ) arrays which span the
!  entire grid

!  set the coordinate differentials, numr, numz and numphi come
!  from the runhydro.h header file
! >>>NOTE<<< the - 3 part comes from getting agreement from scf
!            code and hydrocode.  if numr  = 128 in the scf code
!            that translates to numr = 130 in the hydrocode
dr = 1.0 / (numr - 3.0)
!dr = 1.0 / (numr - 2.0)
dz = dr
dphi = 2.0 * pi * numphiinv
drinv = 1.0 / dr
dzinv = 1.0 / dz
dphiinv = 1.0 / dphi

!  define r array on every processor, use temp here to avoid 
!  coercions from integer to floating types
x = 1.0
offset = 0
do J = rlwb-1,rupb+1
   r(J) = offset + (x - 2.0)*dr
   x = x + 1.0
enddo

!  now define rhf from r
x = 0.5 * dr 
do J = rlwb-1,rupb+1
   rhf(J) = r(J) + x
enddo

!  and define the inverses of both arrays
where(r /= 0.0) 
  rinv = 1.0/r
elsewhere
  rinv = 0.0
endwhere

rhfinv = 1.0/rhf 

! setup the local zhf array
x = 1.0
if( isym /= 1 ) then
   offset = 0
   do K = zlwb-1, zupb+1
      zhf(K) = offset + (x-1.5)*dz
      x = x + 1.0
   enddo
else
   offset = (0*(numz-2) - numz/2)*dz
   do K = zlwb-1,zupb+1
      zhf(K) = offset + (x - 0.5)*dz
      x = x + 1.0
   enddo 
endif

! set up the azimuthal angle
x = 0.0
do L = philwb, phiupb
   phi(L) = x * dphi
   x = x + 1.0
enddo

! global radius array
x = 1.0
do J = 1, numr
   r_g(J) = (x-2.0)*dr
   x = x + 1.0
enddo

! global rhf array
x = 0.5*dr
do J = 1, numr
   rhf_g(J) = r_g(J) + x
enddo

! define the inverse arrays for r and rhf
where( r_g /= 0.0 ) 
   rinv_g = 1.0/r_g
elsewhere
   rinv_g = 0.0
endwhere

rhfinv_g = 1.0/rhf_g

! setup the global zhf array
x = 1.0
if( isym /= 1 ) then
   do K = 1, numz
      zhf_g(K) = (x-1.5)*dz
      x = x + 1.0
   enddo
else
   offset = - (numz/2) * dz
   do K = 1, numz
      zhf_g(K) = offset + (x-0.5)*dz
      x = x + 1.0
   enddo
endif

! trigonemtric arrays for cell centered and vertex centered grid
x = 0.5*dphi
do L = 1, numphi
   cosine(L) = cos(phi(L))
   sine(L) = sin(phi(L))
enddo

call potsetup( have_green_funcs )

return
end subroutine setup
