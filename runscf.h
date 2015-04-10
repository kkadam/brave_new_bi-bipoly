    integer, parameter :: numr = 130
    integer, parameter :: numz = 130
    integer, parameter :: numphi = 256

       real, parameter :: n1 = 1.5
       real, parameter :: nc1 = 1.5
       real, parameter :: mu1 = 1.0
       real, parameter :: muc1 = 1.0

       real, parameter :: n2 = 1.5
       real, parameter :: nc2 = 1.5
       real, parameter :: mu2 = 1.0
       real, parameter :: muc2 = 1.0
       
!       integer, parameter :: phi1 = int(numphi / 4.0) - 1
!       integer, parameter :: phi2 = int(numphi / 4.0) + 1
!       integer, parameter :: phi3 = int(3.0 * numphi / 4.0) - 1
!       integer, parameter :: phi4 = int(3.0 * numphi / 4.0) + 1
       
       integer, parameter :: numr_procs = 1
       integer, parameter :: numz_procs = 1

       integer, parameter :: numr_dd = ( (numr - 2)/numr_procs ) + 2     !numr

       integer, parameter :: numz_dd = ( (numz - 2)/numz_procs ) + 2     !numz 

       integer, parameter :: rlwb = 2, rupb = numr_dd - 1       !2, numr-1

       integer, parameter :: zlwb = 2, zupb = numz_dd - 1       !2, numz-1

       integer, parameter :: philwb = 1, phiupb = numphi

       integer, parameter :: numr_dd_z = (numr-2)/numz_procs + 2

       integer, parameter :: numphi_dd = numphi/numr_procs

       integer, parameter :: numphi_by_two = numphi / 2

       integer, parameter :: numphi_by_four = numphi / 4

       real, parameter :: numphiinv = 1.0 / numphi

       integer, parameter :: maxit = 100 

       real, parameter :: epsilon = 1.0e-5 ! expected minimmum density at the star edge

       real, parameter :: eps = 1e-4  ! Convergance criterion for all parameters

       real, parameter :: densmin = 1e-10 ! Density floor


! restrictions on the above parameters:
!
!  numphi must be a power of two
!
!  (numr - 2) must be evenly divisible by numr_procs
!
!  (numz - 2) must be evenly divisible by numz_procs
!
!  numphi has to be evenly divisible by numr_procs for the
!         data redistribution in the Poisson solver
!
!  (numr-2) must be evenly divisible by numz_procs for the
!           data redistribution in the Poisson solver
!
!  numr_procs must be greater than or equal to two
!
!  numz_procs must be greater than or equal to two
!
!  numr_procs must be an even number
!
!  numz_procs must be an even number
