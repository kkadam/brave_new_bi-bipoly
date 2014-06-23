!************************************************************
!*
!*  MAIN
!*
!************************************************************
program main
  implicit none
  include 'runscf.h'
!************************************************************
!*
!*  Global Variables

  real, dimension(numr,numz,numphi) :: pot, rho
  common /poisson/ pot, rho

  real :: dr, dz, dphi, drinv, dzinv, dphiinv
  common /coord_differentials/ dr, dz, dphi,                   &    
                             drinv, dzinv, dphiinv

  integer :: isym
  integer, dimension(3) :: boundary_condition
  common /boundary_conditions/ isym, boundary_condition

!*
!************************************************************      
!*
!*   Local variables

  logical :: have_green_funcs
  integer :: one, two
  integer :: I, J, N, frnum
  real :: time1, time2
  integer :: qfinal
  integer :: first_model, last_model
  integer :: ierror
  integer :: ra, za, phia
  integer :: rb, zb, phib
  integer :: rc, zc, phic
  integer :: rd, re
  integer :: model_number
  integer :: initial_model_type
  real :: rhom1, rhom2
  real :: frac
  real :: eps

!*
!************************************************************      

! TIMING
  call cpu_time(time1)

  have_green_funcs = .false.
  call setup( have_green_funcs )

  call cpu_time(time2)

  write(*,*) ' setup done in time: ', time2 - time1

  open(unit=20, file='init', form='formatted', status='old')
  read(20,*) first_model, last_model

  za = zlwb
  zb = zlwb
  zc = zlwb
  phia = 1
  phib = 1
  phic = numphi / 2 + 1

  eps = 1.0e-5

  qfinal = 1

  do I = first_model, last_model

    read(20,*) model_number, ra, rb, rc ,rd ,re , initial_model_type, rhom1, rhom2, frac

    call binary_initialize(rhom1, rhom2, ra, rb, rc, phia, phib, phic, &
                            za,  zb, zc, initial_model_type)

    call binary_scf(model_number, initial_model_type, ra, rb, rc, rd, re, rhom1, rhom2, frac, &
                   eps,  qfinal)

    write(*,*) 'Finished Model: ', I
    write(*,*) model_number, ra, rb, rc, rd, re, initial_model_type, rhom1, rhom2, frac

  enddo

  close(20)

  stop
end program main
