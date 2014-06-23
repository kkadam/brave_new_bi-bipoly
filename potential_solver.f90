!*************************************************************
!*
!*  POTENTIAL_SOLVER
!*
!*************************************************************
subroutine potential_solver
implicit none
include 'runscf.h'
!*************************************************************
!*
!  potential_solver is the driver porgram for solving
!  Poisson's equation
!* 
!*************************************************************
!*
!*   Global Variables

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr,numz,numphi) :: potp, rhop
common /potarrays/ potp, rhop

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

!*
!*************************************************************
!*
!*  Local Variables

integer :: nsteps

integer :: J, K, L

!*
!*************************************************************
!  initialize the local variable
nsteps = 0

! copy the density array into a workspace copy
do L = 1, numphi
   do K = 1, numz
      do J = 1, numr
         rhop(J,K,L) = rho(J,K,L)
      enddo
   enddo
enddo
! need to zero out sections of rhop to avoid problems
! with material piling up in the boundary zones for
! the dirichlet boundary conditions

! Top of the grid
do L = 1, numphi
   do K = zupb, zupb+1
      do J = 1, numr
         rhop(J,K,L) = 0.0
      enddo
   enddo
enddo
! outer edge of grid
do L = 1, numphi
   do K = 1, numz
      do J = rupb, rupb+1
         rhop(J,K,L) = 0.0
      enddo
   enddo
enddo
! bottom of grid
if( isym == 1 ) then
   do L = 1, numphi
      do K = zlwb-1, zlwb
         do J = 1, numr
            rhop(J,K,L) = 0.0
         enddo
      enddo
   enddo
endif

if( tstep > 1 ) then
   ! have a previous set of values for the potential
   ! it has either been read in from the continuation
   ! file by setup or exists from the end of the last
   ! time potential_solver was called.  Perform only
   ! 5 ADI iterations as we are starting from a good
   ! guess for the potential
   nsteps = 20
   !nsteps = 10
else
   ! don't have a guess for the potential but 20 ADI
   ! iterations are sufficient for a cold start to
   ! get an acceptable solution
   !nsteps = 50
   nsteps = 50
   do L = 1, numphi
      do K = 1, numz
         do J = 1, numr
            potp(J,K,L) = 0.0
         enddo
      enddo
   enddo
endif

! solve for the boundary values of the potential
call bessel 

! solve Poisson's equation
call helmadi(nsteps)
  
! fill in the potential with the solution
! if you want to add an external potential in
! addition to the self gravity this would be a
! good place to do it
do L = 1, numphi
   do K = 1, numz
      do J = 1, numr
         pot(J,K,L) = potp(J,K,L)
      enddo
   enddo
enddo

return
end subroutine potential_solver
