!************************************************************
!*
!*  MAIN
!*
!************************************************************
program main
implicit none
include 'runscf.h'
!include 'mpif.h'
!************************************************************
!*
!*  Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi,                   &    
                             drinv, dzinv, dphiinv

integer, dimension(3) :: boundary_condition
common /boundary_conditions/ boundary_condition

logical :: iam_on_top, iam_on_bottom, iam_on_axis,           &
           iam_on_edge, iam_root
integer :: column_num, row_num
integer :: iam, down_neighbor, up_neighbor,                  &
           in_neighbor, out_neighbor, root,                  &
           REAL_SIZE, INT_SIZE, numprocs
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,           &
                        iam_on_bottom, iam_on_axis,          &
                        iam_on_edge, down_neighbor,          &
                        up_neighbor, in_neighbor,            &
                        out_neighbor, root, column_num,      &
                        row_num, pe_grid, iam_root,          &
                        REAL_SIZE, INT_SIZE

!*
!************************************************************      
!*
!*   Local variables

integer, dimension(numr_procs*numz_procs,2) :: pe_coord

integer :: one, two

integer :: I, J, N

real :: time1, time2

integer :: qfinal

integer :: first_model, last_model

integer :: ra, za, phia

integer :: rb, zb, phib

integer :: rc, zc, phic

integer :: rd, re

integer :: model_number

integer :: initial_model_type

real :: rhom1, rhom2

real :: frac

logical ::have_green_funcs
!*
!************************************************************      

! set the default data sizes for MPI communications
REAL_SIZE = 8
INT_SIZE = 8

! TIMING
call cpu_time(time1)

! Initialize local variables
pe_coord = 0

!  Set up logical and integer variables that describe the
!  processor grid and the message passing pattern

iam_on_top = .false.
iam_on_bottom = .false.
iam_on_axis = .false.
iam_on_edge = .false.
iam_root = .false.
root = 0

iam = 0
numprocs = numr_procs * numz_procs
if( iam == 0 ) iam_root = .true.

!  Make sure the input data in runhdro.h doesn't violate
!  assumptions made in the code implementaion.  If there
!  is a problem might as well stop now...
!if( mod(numr_procs,2) /= 0 ) then
!   ! not en even number of processors in radial
!   ! direction to avoid deadlock in message passing
!   if( iam_root ) write(6,*) 'Stoppinng, bad numr_procs'
!   stop
!endif

!if( mod(numz_procs,2) /= 0 ) then
!   ! not an even number of processors in vertical
!   ! direction to avoid deadlock in message passing
!   if( iam_root ) write(6,*) 'Stopping, bad numz_procs'
!   stop
!endif

!if( mod(numphi,numr_procs) /= 0 ) then
!   ! to do data swapping in helmadi the azimuthal
!   ! dimension of the data has to distribute evenly
!   ! across numr_procs
!   if( iam_root ) write(6,*) 'Stopping, wrong numphi, numr_procs'
!   stop
!endif

!if( mod(numr-2,numz_procs) /= 0 ) then
!   ! helmadi's communications also require that
!   ! radial data be divisible
!   ! evenly by numz_procs.  Checking if I picked
!   ! a valid pair of numr and numz_procs
!   if( iam_root ) write(6,*) 'Stopping, wrong numr, numz_procs'
!   stop
!endif

!  Setup the logical pe grid
N = 1
do J = 1, numz_procs
   do I = 1, numr_procs
      pe_coord(N,1) = I-1
      pe_coord(N,2) = J-1
      pe_grid(I,J) = N-1
      N = N + 1
   enddo
enddo

!  Determine which processors have external boundary zones
!  in place of a guard cell layer that will require special
!  treatment
if( pe_coord(iam+1,2) == 0 ) iam_on_bottom = .true.

if( pe_coord(iam+1,2) == numz_procs-1 ) iam_on_top = .true.

if( pe_coord(iam+1,1) == 0 ) iam_on_axis = .true.

if( pe_coord(iam+1,1) == numr_procs-1 ) iam_on_edge = .true. 

!  Determine passing partners for comm and related routines
if( iam_on_top ) then
   up_neighbor = -1
else
   up_neighbor = iam + numr_procs
endif

if( iam_on_bottom ) then
   down_neighbor = -1
else
   down_neighbor = iam - numr_procs
endif

if( iam_on_axis ) then
   in_neighbor = -1
else
   in_neighbor = iam - 1
endif

if( iam_on_edge ) then
   out_neighbor = -1
else
   out_neighbor = iam + 1
endif

!  Determine whether pe is in an even row or column to order
!  message passing
row_num = pe_coord(iam+1,2)
 
column_num = pe_coord(iam+1,1)


 have_green_funcs = .true.
call setup( have_green_funcs )
 have_green_funcs = .true.
call cpu_time(time2)
 have_green_funcs = .true.
if ( iam_root ) then
   write(*,*) ' setup done in time: ', time2 - time1
endif


open(unit=20, file='init', form='formatted', status='old')
read(20,*) first_model, last_model

za = zlwb
zb = zlwb
zc = zlwb
phia = 1
phib = 1
phic = numphi / 2 + 1

qfinal = 1

do I = first_model, last_model

   read(20,*) model_number, ra, rb, rc ,rd ,re , initial_model_type, rhom1, rhom2, frac

   call binary_initialize(rhom1, rhom2, ra, rb, rc, phia, phib, phic, &
                            za,  zb, zc, initial_model_type)
print*, "Finished initialize"
   call binary_scf(model_number, initial_model_type, ra, rb, rc, rd, re, rhom1, rhom2, frac, qfinal)

   write(*,*) 'Finished Model: ', I
   write(*,*) model_number, ra, rb, rc, initial_model_type, rhom1, rhom2, frac

enddo

close(20)


stop
end program main
