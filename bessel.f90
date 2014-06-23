!****************************************************************
!*
!*   BESSEL
!*
!****************************************************************
subroutine bessel
implicit none
include 'runscf.h'
include 'pot.h'
!****************************************************************
!*
!  bessel calculates the boundary values of the gravitational
!  potential at
! 
!  -> top surface, K = numz_dd for all J and L for pes on
!                  top of the grid
!
!  -> side surface, J = numr_dd for all K and L for pes on
!                   side of the grid
!
!  -> bottom surface if isym = 1 (no assumed symmetry), K =1
!                    for all J and L for pes on bottom of
!                    the grid
!    isym = 1 is not supported.
!
!  Boundary potential values are calculated by convolving the
!  density distribution with the appropriate cylindrical Green
!  function.  Initially implemented by Howard Cohl.  See his
!  thesis for discussion and orignal hpf source code.
!
!*
!****************************************************************
!*
!*  Global Variables

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr,numz,numphi) :: potp, rhop
common /potarrays/ potp, rhop

real, dimension(numr,numz,numr,mmax) :: tmr
real, dimension(numr,numz,numz,mmax) :: smz
common /green_functions/ tmr, smz

real, dimension(numphi,mmax) :: bes_cos, bes_sin
common /bessel_trig/ bes_cos, bes_sin

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi,                             &
                             drinv, dzinv, dphiinv
    
integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

!*
!****************************************************************
!*
!*  Local Variables

real, dimension(numr,numz) :: TMPC, TMPS

real, dimension(numr,numphi) :: phitTMP, pott, phitTMPC,            &
                                   phitTMPS

real, dimension(numz,numphi) :: phisTMP, pots, phisTMPC,            &
                                   phisTMPS

real, dimension(numr,mmax) :: StC, StS

real, dimension(numz,mmax) :: SsC, SsS

real :: factor

integer :: J, K, L, M, lwrb, uprb, counter, rindex, zindex

integer :: I

!*
!****************************************************************
! initialize the local variables
do K = 1, numz
   do J = 1, numr
      TMPC(J,K) = 0.0
      TMPS(J,K) = 0.0
   enddo
enddo
do L = 1, numphi
   do J = 1, numr
      phitTMP(J,L)  = 0.0
      pott(J,L)     = 0.0
      phitTMPC(J,L) = 0.0
      phitTMPS(J,L) = 0.0
   enddo
enddo
do L = 1, numphi
   do K = 1, numz
      phisTMP(K,L)  = 0.0
      pots(K,L)     = 0.0
      phisTMPC(K,L) = 0.0
      phisTMPS(K,L) = 0.0
   enddo
enddo
do M = 1, mmax
   do J = 1, numr
      StC(J,M) = 0.0
      StS(J,M) = 0.0
   enddo
enddo
do M = 1, mmax
   do K = 1, numz
      SsC(K,M) = 0.0
      SsS(K,M) = 0.0
   enddo
enddo
lwrb = 0
uprb = 0
counter = 0

!  factor is the common multiplier for converting summation of
!  Green function times the density to a potential
if( isym == 3 ) then
   factor = - 2.0 * dr * dz * dphi
else
   factor = - dr * dz * dphi
endif

!  evaluate the m=0 contribution to top and side slices
!  of the potential as a special case because the sine terms 
!  drop out.
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         TMPC(J,K) = TMPC(J,K) + rho(J,K,L)
      enddo
   enddo
enddo
do J = 2, numr-1
   do zindex = zlwb, zupb
      do rindex = rlwb, rupb
         StC(J,1) = StC(J,1) + tmr(rindex,zindex,J,1) *             &
                               TMPC(rindex,zindex) 
      enddo
   enddo
enddo 
do K = 2, numz-1 
   do zindex = zlwb, zupb
      do rindex = rlwb, rupb
         SsC(K,1) = SsC(K,1) + smz(rindex,zindex,K,1) *             &
                               TMPC(rindex,zindex)
      enddo
   enddo 
enddo

!  now compute the contributions to the boundary potential for
!  modes with m > 0 up to mmax - 1
do M = 2, mmax
   TMPC = 0.0
   TMPS = 0.0
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            TMPC(J,K) = TMPC(J,K) + rho(J,K,L)*bes_cos(L,M)
            TMPS(J,K) = TMPS(J,K) + rho(J,K,L)*bes_sin(L,M)
         enddo
      enddo
   enddo
   do J = 2, numr-1
      do zindex = zlwb, zupb
         do rindex = rlwb, rupb
            StC(J,M) = StC(J,M) + tmr(rindex,zindex,J,M) *           &
                                  TMPC(rindex,zindex)
            StS(J,M) = StS(J,M) + tmr(rindex,zindex,J,M) *           &
                                  TMPS(rindex,zindex)
         enddo
      enddo
   enddo
   do K = 2, numz-1
      do zindex = zlwb, zupb
         do rindex = rlwb, rupb
            SsC(K,M) = SsC(K,M) + smz(rindex,zindex,K,M) *           &
                                  TMPC(rindex,zindex)
            SsS(K,M) = SsS(K,M) + smz(rindex,zindex,K,M) *           &
                                  TMPS(rindex,zindex)
         enddo
      enddo
   enddo
enddo

! reduce the convolution
! of G(r|r') with rho to a potential at the
! top of the grid
do L = philwb, phiupb
   do J = rlwb, rupb
      phitTMP(J,L) = StC(J,1)
   enddo
enddo
do M = 2, mmax
   do L = philwb, phiupb
      do J = rlwb, rupb
         phitTMPC(J,L) = StC(J,M)*bes_cos(L,M)
         phitTMPS(J,L) = StS(J,M)*bes_sin(L,M)
      enddo
   enddo
   do J = rlwb, rupb
      phitTMP(J,:) = phitTMP(J,:) + 2.0*phitTMPC(J,:) +             &
                     2.0*phitTMPS(J,:)
   enddo
enddo
do L = philwb, phiupb
   do J = rlwb, rupb
      pott(J,L) = factor * phitTMP(J,L)
   enddo
enddo
if( isym == 3 ) then
   do L = philwb, phiupb
      pott(1,L) = pott(2,L)
   enddo
else
   do L = 1, numphi_by_two
      pott(1,L) = pott(2,L+numphi_by_two)
      pott(1,L+numphi_by_two) = pott(2,L)
   enddo
endif
do L = 1, numphi
   do J = 1, numr
      potp(J,numz,L) = pott(J,L)
   enddo
enddo
! done calculating top boundary potential

! reduce the convolution
! of G(r|r') with rho to a potential at the
! outer edge of the grid
do L = philwb, phiupb
   do K = zlwb, zupb
      phisTMP(K,L) = SsC(K,1)
   enddo
enddo
do M = 2, mmax
   do L = philwb, phiupb
      do K = zlwb, zupb
         phisTMPC(K,L) = SsC(K,M)*bes_cos(L,M)
         phisTMPS(K,L) = SsS(K,M)*bes_sin(L,M)
      enddo
   enddo
   do L = philwb, phiupb
      do K = zlwb, zupb
         phisTMP(K,L) = phisTMP(K,L) + 2.0*phisTMPC(K,L) +              &
                        2.0*phisTMPS(K,L)
      enddo
   enddo
enddo
do L = philwb, phiupb
   do K = zlwb, zupb
      pots(K,L) = factor*phisTMP(K,L)
   enddo
enddo
if( isym /= 1 )  then
   do L = 1, numphi
      pots(1,L) = pots(2,L)
   enddo
endif
do L = 1, numphi
   do K = 1, numz
      potp(numr,K,L) = pots(K,L)
   enddo
enddo
! done calculating side boundary potential

return
end subroutine bessel
