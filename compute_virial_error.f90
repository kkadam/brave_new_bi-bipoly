subroutine compute_virial_error(psi, h, omega, pin, volume_factor, virial_error1, virial_error2, virial_error)
implicit none
include 'runscf.h'
!include 'mpif.h'
!*****************************************************************************************
!
!  subroutine arguments
!

real, dimension(numr_dd, numphi), intent(in) :: psi

real, dimension(numr_dd,numz_dd,numphi) :: h

real, intent(in) :: omega

real, intent(in) :: pin

real, intent(in) :: volume_factor

real, intent(out) :: virial_error1

real, intent(out) :: virial_error2

real, intent(out) :: virial_error

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

real, dimension(numr_dd,numz_dd,numphi) :: temp

real :: w1,  w2, wtot

real :: t1, t2, ttot

real :: s1, s2, stot

real :: ret1, ret2, global_ret1, global_ret2

real :: gamma

integer :: phi1, phi2, phi3, phi4

integer :: I, J, K

integer :: ierror

!
!*****************************************************************************************

virial_error = 0.0
virial_error1 = 0.0
virial_error2 = 0.0

phi1 = int(numphi / 4.0) - 1
phi2 = int(numphi /  4.0) + 1
phi3 = int(3.0 * numphi / 4.0) - 1
phi4 = int(3.0 * numphi / 4.0) + 1

gamma = 1.0 + 1.0 / pin

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

virial_error  = abs(2.0*ttot + 3.0*stot + wtot) / abs(wtot)
virial_error1 = abs(2.0*t1   + 3.0*s1   + w1  ) / abs(w1)
virial_error2 = abs(2.0*t2   + 3.0*s2   + w2  ) / abs(w2)

end subroutine compute_virial_error
