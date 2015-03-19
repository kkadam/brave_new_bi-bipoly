subroutine ancient_output(c1, c2, omsq, hm1, hm2, mass1, mass2, psi, h,           &
                         qfinal, initial_model_type, model_number, ra, za, phia, &
                         rb, zb,  phib, rc, zc, phic, rhm1, rhm2, pin, rhom1,    &
                         rhom2, xavg1, xavg2, separation,  com, volume_factor,   &
                         eps)
implicit none
include 'runscf.h'
!include 'mpif.h'
!*****************************************************************************************
!
!  subroutine arguments
!

real, dimension(maxit), intent(in) :: c1, c2, omsq, hm1, hm2, mass1, mass2

real, dimension(numr_dd, numphi), intent(in) :: psi

real, dimension(numr_dd,numz_dd,numphi) :: h

integer, intent(in) :: qfinal

integer, intent(in) :: initial_model_type

integer, intent(in) :: ra, za, phia

integer, intent(in) :: rb, zb, phib

integer, intent(in) :: rc, zc, phic

integer, intent(in) :: model_number

real, dimension(3), intent(in) :: rhm1, rhm2

real, intent(in) :: pin

real, intent(in) :: rhom1, rhom2

real, intent(in) :: xavg1, xavg2, separation, com

real, intent(in) :: volume_factor

real, intent(in) :: eps

!
!*****************************************************************************************
!
!  global variables
!

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numphi) :: cosine, sine
common /trig/ cosine, sine

real :: pi, grav
common /constants/ pi, grav

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

!
!*****************************************************************************************
!
! locall variables
!

real, dimension(numr_dd,numz_dd,numphi) :: rchpot

real, dimension(numr_dd,numz_dd,numphi) :: temp

real, parameter :: epsilon = 1.0e-5 ! expected minimmum density at the star edge

real :: gamma

real :: volr1, volr2, reffr1, reffr2, global_volr1, global_volr2

real :: vol1, vol2, reff1, reff2

real :: yavg1, yavg2

real :: en1, en2, entot

real :: e1, e2, etot

real :: j1, j2, jtot

real :: w1,  w2, wtot

real :: t1, t2, ttot

real :: s1, s2, stot

real :: virialerr1, virialerr2, virialerr

real :: kappa1, kappa2

real :: pm1, pm2

real :: period, omega

real :: kepler

real :: ret1, ret2, global_ret1, global_ret2

real :: rpotcrit, xcrit, cuurvature, rchtest, rchmax,  rchmin

real :: temp_rpotcrit, temp_xcrit, curvature

integer :: isave, flag

real :: my_rchmin, my_rchmax

real, dimension(3) :: rmaxloc, rminloc, my_rmaxloc, my_rminloc

integer :: primary

real :: star2maxr, temp_star2maxr

real :: rochemax1, rochemax2, temp_rochemax1, temp_rochemax2

real :: temp_rch

real, dimension(3) :: temp_rch_loc

real :: l2loc, l3loc, temp_l2loc, temp_l3loc

integer :: louter1,  louter2

integer :: phi1, phi2, phi3, phi4

integer :: I, J, K, L

integer :: index

integer :: ierror

!integer, dimension(MPI_STATUS_SIZE) :: istatus

character(len=50) :: model_template

character(len=56) :: model_file

!
!*****************************************************************************************

model_template = 'ancient_details_'

phi1 = int(numphi / 4.0) - 1
phi2 = int(numphi /  4.0) + 1
phi3 = int(3.0 * numphi / 4.0) - 1
phi4 = int(3.0 * numphi / 4.0) + 1

gamma = 1.0 + 1.0 / pin

primary = 1
if ( mass2(qfinal) > mass1(qfinal) ) then
   primary = 2
endif

omega = sqrt(omsq(qfinal))

! sum up the virial pressuure
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * h(I,J,K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
!call mpi_reduce(ret1, global_ret1, 1,  REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1,  REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1,  REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1,  REAL_SIZE, root, MPI_COMM_WORLD, ierror)
s1 = volume_factor * ret1 / (pin+1.0)
s2 = volume_factor * ret2 / (pin+1.0)
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
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
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
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
t1 = - omega * omega * volume_factor * ret1
t2 = - omega * omega * volume_factor * ret2
ttot = t1 + t2

virialerr  = abs(2.0*ttot + 3.0*stot + wtot) / abs(wtot)
virialerr1 = abs(2.0*t1   + 3.0*s1   + w1  ) / abs(w1)
virialerr2 = abs(2.0*t2   + 3.0*s2   + w2  ) / abs(w2)

pm1 = rhom1 * hm1(qfinal) / (pin + 1.0)
pm2 = rhom2 * hm2(qfinal) / (pin + 1.0)
kappa1 = pm1 / rhom1**gamma
kappa2 = pm2 / rhom2**gamma

period = 2.0 * pi /  omega
kepler = (separation**3) * omega * omega / (mass1(qfinal) + mass2(qfinal))

! sum the angular momentum
do K = philwb, phiupb
   do  J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * psi(I,K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
j1 = - 2.0 * omega * volume_factor * ret1
j2 = - 2.0 * omega * volume_factor * ret2
jtot = j1 + j2

! sum the internal energiies
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * rho(I,J,K)**gamma
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
e1 = pin * volume_factor * kappa1 * ret1
e2 = pin * volume_factor * kappa2 * ret2
etot = e1 + e2

! total ennergies
en1 = t1 + e1 + w1
en2 = t2 + e2 + w2
entot = ttot + etot + wtot

! sum the total volume for each star
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         if ( rho(I,J,K) > 0.0 ) then
            temp(I,J,K) = rhf(I)
         else
            temp(I,J,K)  = 0.0
         endif
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
vol1 = volume_factor * ret1
vol2 = volume_factor * ret2
reff1 = (0.75 * vol1 / pi)**(1.0/3.0)
reff2 = (0.75 * vol2 / pi)**(1.0/3.0)

! calculate the y moment oof the density distribution
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * rhf(I) * sine(K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)
!call mpi_reduce(ret1, global_ret1, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_reduce(ret2, global_ret2, 1, REAL_SIZE, MPI_SUM, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret1, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
!call mpi_bcast(global_ret2, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
yavg1 = volume_factor * ret1 / mass1(qfinal)
yavg2 = volume_factor * ret2 / mass2(qfinal)

! compute the Roche potential
do K = philwb, phiupb
   do J = zlwb-1, zupb+1
      do I = rlwb-1, rupb+1
         rchpot(I,J,K) = pot(I,J,K) + omega * omega * psi(I,K)
      enddo
   enddo
enddo

  do L = 1, numphi_by_two
    do K = 1, numz
      rchpot(rlwb-1,K,L)               = rchpot(rlwb,K,L+numphi_by_two) !wth?
      rchpot(rlwb-1,K,L+numphi_by_two) = rchpot(rlwb,K,L)!rho(rlwb,K,L)
    enddo
  enddo
! find  the minimum value of the Roche potential  wth, can be fixed
  rchmin = 0.0
  do K = philwb, phiupb
    do J = zlwb, zupb
      do I = rlwb, rupb
        if ( rchpot(I,J,K) < rchmin ) then
           rchmin = rchpot(I,J,K)
           rminloc(1) = rhf(I)
           rminloc(2) = zhf(J)
           rminloc(3) = phi(K)
        endif
      enddo
    enddo
  enddo

! find the maximum value of the Roche potential
  rchmax = - 1.0e6
  do K = philwb, phiupb
    do J = zlwb, zupb
      do I = rlwb, rupb
        if ( rchpot(I,J,K) > rchmax ) then
          rchmax = rchpot(I,J,K)
          rmaxloc(1) = rhf(I)
          rmaxloc(2) = zhf(J)
          rmaxloc(3) = phi(K)
        endif
      enddo
    enddo
  enddo

! find the location of the L1 point and the potenttial there
flag = 0
isave = 0
xcrit = 1.0
rpotcrit = 1.0
do I = rupb, rlwb, -1
   if ( rhf(I) < rhm2(1) ) then
      rchtest = (rchpot(I,zlwb,phic) - rchpot(I+1,zlwb,phic)) * &
                (rchpot(I-1,zlwb,phic) - rchpot(I,zlwb,phic))
      if ( rchtest < 0.0 ) then
         curvature = rchpot(I+1,zlwb,phic) + rchpot(I-1,zlwb,phic) - 2.0 * rchpot(I,zlwb,phic)
         if ( cuurvature < 0.0 ) then
            xcrit = - rhf(I)
            rpotcrit = rchpot(I,zlwb,phic)
            isave = I
            flag = 0
         endif
      endif
   endif
enddo
! L1 is not on the -ve x axis if isave is zero
if ( isave == 0 ) then
   do I = rlwb, rupb
       if  ( rhf(I) < rhm1(1) ) then
          rchtest = (rchpot(I+1,zlwb,phia) - rchpot(I,zlwb,phia)) * &
                    (rchpot(I,zlwb,phia) - rchpot(I-1,zlwb,phia))
          if ( rchtest < 0.0 ) then
             curvature = rchpot(I+1,zlwb,phia) + rchpot(I-1,zlwb,phia) - 2.0 * rchpot(I,zlwb,phia)
             if ( curvature < 0.0 ) then
               xcrit = rhf(I)
               rpotcrit = rchpot(I,zlwb,phia)
               isave = I
               flag = 1
            endif
          endif
      endif
   enddo
endif
if ( isave == 0 ) then
   ! the  L1 point  muust bbe on the axxis
   xcrit = 0.0
   rpotcrit = 0.5 * ( rchpot(rlwb-1,zlwb,phia) + rchpot(rlwb,zlwb,phia) )
endif

! find the L2 and  L3 points if they are on the computational grid
l2loc = 0.0
l3loc = 0.0
do I = rlwb, rupb
   if ( rhf(I) >  rhm1(1) ) then
       rchtest = (rchpot(I+1,zlwb,phia) - rchpot(I,zlwb,phia)) * &
                 (rchpot(I,zlwb,phia) - rchpot(I-1,zlwb,phia))
       if ( rchtest < 0.0 ) then
          curvature = rchpot(I+1,zlwb,phia) + rchpot(I-1,zlwb,phia) - 2.0  * rchpot(I,zlwb,phia)
          if ( curvature < 0.0 ) then
              l2loc = rhf(I)
              exit
          endif
       endif
   endif
enddo
do I = rlwb, rupb
   if ( rhf(I) > rhm2(1) ) then
      rchtest = (rchpot(I+1,zlwb,phic) - rchpot(I,zlwb,phic)) * &
                (rchpot(I,zlwb,phic) - rchpot(I-1,zlwb,phic))
      if ( rchtest < 0.0 ) then
         curvature = rchpot(I+1,zlwb,phic) + rchpot(I-1,zlwb,phic) - 2.0  * rchpot(I,zlwb,phic)
         if ( curvature < 0.0 ) then
            l3loc =  rhf(I)
            exit
         endif
      endif
   endif
enddo

! find the outer edge of the Roche lobes
do I = rlwb, rupb
   if ( rhf(I) > rhm1(1) ) then
      if ( rchpot(I,zlwb,phia) <= rpotcrit .and. rchpot(I+1,zlwb,phia) >= rpotcrit ) then
         rochemax1 = rhf(I+1)
         exit
      endif
   endif
enddo
do I = rlwb, rupb
   if ( rhf(I) > rhm2(1) ) then
      if ( rchpot(I,zlwb,phic) <= rpotcrit .and. rchpot(I+1,zlwb,phic) >= rpotcrit ) then
         rochemax2 = rhf(I+1)
         exit
      endif
   endif
enddo

! total up the volume in each Roche lobe
volr1 = 0.0
volr2 = 0.0
if ( xcrit >= 0.0 ) then
   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rchpot(I,J,K) <= rpotcrit ) then
               if ( rhf(I) * cosine(K) < xcrit ) then
                  volr2 = volr2 + rhf(I)
               else
                  volr1 = volr1 + rhf(I)
               endif
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rchpot(I,J,K) <= rpotcrit ) then
               if ( rhf(I) * cosine(K) < xcrit ) then
                  volr2 = volr2 + rhf(I)
               else
                  volr1 = volr1 + rhf(I)
               endif
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rchpot(I,J,K) <= rpotcrit ) then
               volr2 = volr2 + rhf(I)
            endif
         enddo
      enddo
   enddo
else
   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rchpot(I,J,K) <= rpotcrit ) then
               volr1 = volr1 + rhf(I)
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rchpot(I,J,K) <= rpotcrit ) then
               volr1 = volr1 + rhf(I)
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rchpot(I,J,K) <= rpotcrit ) then
               if ( rhf(I) * cosine(K) > xcrit ) then
                  volr1 = volr1 + rhf(I)
               else
                  volr2 = volr2 + rhf(I)
               endif
            endif
         enddo
      enddo
   enddo
endif
volr1 = volume_factor * volr1
volr2 = volume_factor * volr2
reffr1 = (0.75 * volr1 / pi)**(1.0/3.0)
reffr2 = (0.75 * volr2 / pi)**(1.0/3.0)

! find the outer edge of the star on the -ve x axis
star2maxr = 0.0
do I = rlwb, rupb
   if  ( rho(I,zlwb,phic) > epsilon .and. rho(I+1,zlwb,phic) < epsilon ) then
      star2maxr = rhf(I)
   endif
enddo

if ( iam_root ) then

   write(model_file,'(a,i6)') trim(model_template), model_number
   open(unit=11, file=trim(model_file),form='formatted',status='unknown')

   write(11,*)  'Model Number: ', model_number
   write(11,*) 
   write(11,*) 'For Star 1:'
   write(11,*) 'Mass 1: ', mass1(qfinal)
   write(11,*) '(<x>,<y>): ', xavg1, yavg1
   write(11,*) 'Maximum Density: ', rhom1
   write(11,*) 'Polytropic Constant: ', kappa1
   write(11,*) 'Virial Pressure: ', s1
   write(11,*) 'Potential Energy: ', w1
   write(11,*) 'Kiinetic Energy: ', t1
   write(11,*) 'Virial Error: ', virialerr1
   write(11,*) 'Pressure Maximum: ', pm1
   write(11,*) 'Enthaply Maximum: ', hm1(qfinal)
   write(11,*) 'Maximum at (r, z, phi): ', rhm1
   write(11,*) 'Inner Boundary Point: ', rb, zb, phib
   write(11,*) '(r, z, phi): ', rhf_g(rb), zhf_g(zb), phi(phib)
   write(11,*) 'Outer Boundary Point: ', ra, za, phia
   write(11,*) '(r, z, phi): ', rhf_g(ra), zhf_g(za), phi(phia)
   write(11,*) 'Volume: ', vol1
   write(11,*) 'Effective Radius: ', reff1
   write(11,*) 'Roche Volume: ', volr1
   write(11,*)  'Effective Roche Radius: ', reffr1
   write(11,*) 'Angular Momentum: ', j1
   write(11,*) 'Internal Energy: ', e1
   write(11,*) 'Total Energy: ', en1
   write(11,*) 'Outer Lagrange Point: ', l2loc
   write(11,*) 'Outer Edge of Roche Lobe: ', rochemax1
   write(11,*)
   write(11,*) 'For Star 2:'
   write(11,*) 'Mass 2: ', mass2(qfinal)
   write(11,*) '(<x>,<y>): ', xavg2, yavg2
   write(11,*) 'Maximum Density: ', rhom2
   write(11,*) 'Polytropic Constant: ', kappa2
   write(11,*) 'Virial Pressure: ', s2
   write(11,*) 'Potential Energy: ', w2
   write(11,*) 'Kinetic Energy: ', t2
   write(11,*) 'Virial Error: ', virialerr2
   write(11,*) 'Pressure Maximum: ', pm2
   write(11,*) 'Enthalpy Maximum: ', hm2(qfinal)
   write(11,*) 'Maximum at (r, z, phi): ', rhm2
   write(11,*) 'Inner Boundary Point: ', rc, zc, phic
   write(11,*) '(r, z, phi): ', rhf_g(rc), zhf_g(zc), phi(phic)
   write(11,*) 'Volume: ', vol2
   write(11,*) 'Effective Radius: ', reff2
   write(11,*) 'Roche Volume: ', volr2
   write(11,*) 'Efffective Roche Radius: ', reffr2
   write(11,*) 'Star outer extent: ', star2maxr
   write(11,*) 'Angular Momentum: ', j2
   write(11,*) 'Internal Energy: ', e2
   write(11,*) 'Total Energy: ', en2
   write(11,*) 'Outer Lagrannge Point: ', l3loc
   write(11,*) 'Outer Edge of Roche Lobe: ', rochemax2
   write(11,*)
   write(11,*) 'Mass Ratio: ', mass1(qfinal) / mass2(qfinal)
   write(11,*) 'Primary is: ', primary
   write(11,*) 'Virial Error: ', virialerr
   write(11,*) 'Center of Mass: ', com
   write(11,*) 'Separation: ', separation
   write(11,*) 'Angular Frequency: ', omega
   write(11,*) 'Period: ', period
   write(11,*) 'Keplers 3rd Constant: ', kepler
   write(11,*) 'Integration Constant 1: ', c1(qfinal)
   write(11,*) 'Integration Constant 2: ', c2(qfinal)
   write(11,*) 'Total Angular Momentum: ', jtot
   write(11,*) 'Total Energy: ', entot
   write(11,*) 'Roche Potential at  L1: ', rpotcrit
   write(11,*) 'x Coordinate of L1: ', xcrit
   write(11,*) 'Maximum Value of Roche Potential: ', rchmax
   write(11,*) 'Located at: ', rmaxloc
   write(11,*) 'Minimum Value of Rohhe Potential: ', rchmin
   write(11,*) 'Located at: ', rminloc
   write(11,*) 'Convergence Criterion: ', eps
   write(11,*) 'Number of Iterations: ', qfinal-1
   write(11,*) 'Polytropic Index: ', pin
   write(11,*) 'Polytropic Exponent: ', gamma
   write(11,*) 'Initial Model Type: ', initial_model_type
   write(11,*)
   close(11)

endif

end subroutine ancient_output
