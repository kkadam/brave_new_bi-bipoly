!****************************************************************
!*
!*   BESSEL
!*
!****************************************************************
subroutine bessel
implicit none
!include 'mpif.h'
include 'runscf.h'
include 'pot.h'
!****************************************************************
!*
!  bessel calculates the boundary values of the gravitational
!  potential at
! 
!  -> top surface, K = numz for all J and L for pes on
!                  top of the grid
!
!  -> side surface, J = numr for all K and L for pes on
!                   side of the grid
!
!  -> bottom surface if isym = 1 (no assumed symmetry), K =1
!                    for all J and L for pes on bottom of
!                    the grid
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
    
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ boundary_condition

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

real, dimension(numr,mmax) :: sum_top_C, sum_top_S

real, dimension(numz,mmax) :: sum_sid_C, sum_sid_S

real, dimension(numr-2,mmax) :: t_buff_C, t_buff_S

real, dimension(numr-2,mmax) :: t_buff_C_summed, t_buff_S_summed

real, dimension(numz-2,mmax) :: s_buff_C, s_buff_S

real, dimension(numz-2,mmax) :: s_buff_C_summed, s_buff_S_summed

real :: factor

integer :: J, K, L, M, lwrb, uprb, counter, rindex, zindex

integer :: I

!*
!****************************************************************
! initialize the local variables
TMPC = 0.0
TMPS = 0.0
phitTMP = 0.0
pott = 0.0
phitTMPC = 0.0
phitTMPS = 0.0
phisTMP = 0.0
pots = 0.0
phisTMPC = 0.0
phisTMPS = 0.0
StC = 0.0
StS = 0.0
SsC = 0.0
SsS = 0.0
sum_top_C = 0.0
sum_top_S = 0.0
sum_sid_C = 0.0
sum_sid_S = 0.0
t_buff_C = 0.0
t_buff_C_summed = 0.0
t_buff_S = 0.0
t_buff_S_summed = 0.0
s_buff_C = 0.0
s_buff_C_summed = 0.0
s_buff_S = 0.0
s_buff_S_summed = 0.0
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

lwrb = 2
uprb = numr - 1

   do M = 1, mmax
      counter = lwrb
      do J = 1, numr - 2
         t_buff_C(J,M) = StC(counter,M)
         t_buff_S(J,M) = StS(counter,M)
         counter = counter + 1
      enddo
   enddo
   t_buff_C_summed = t_buff_C
   t_buff_S_summed = t_buff_S

      do M = 1, mmax
         do J = 1, numr - 2
            sum_top_C(J+1,M) = t_buff_C_summed(J,M)
            sum_top_S(J+1,M) = t_buff_S_summed(J,M)
         enddo
      enddo
   lwrb = uprb + 1
   uprb = lwrb + numr - 3

lwrb = 2
uprb = numz - 1
   do M = 1, mmax
      counter = lwrb
      do K = 1, numz - 2
         s_buff_C(K,M) = SsC(counter,M)
         s_buff_S(K,M) = SsS(counter,M)
         counter = counter + 1
      enddo
   enddo
   s_buff_C_summed = s_buff_C
   s_buff_S_summed = s_buff_S

      do M = 1, mmax
         do K = 1, numz - 2
            sum_sid_C(K+1,M) = s_buff_C_summed(K,M)
            sum_sid_S(K+1,M) = s_buff_S_summed(K,M)
         enddo
      enddo

   lwrb = uprb + 1
   uprb = lwrb + numz - 3

! if on top of the pe grid reduce the convolution
! of G(r|r') with rho to a potential at the
! top of the grid
   do L = philwb, phiupb
      do J = rlwb, rupb
         phitTMP(J,L) = sum_top_C(J,1)
      enddo
   enddo
   do M = 2, mmax
      do L = philwb, phiupb
         do J = rlwb, rupb
            phitTMPC(J,L) = sum_top_C(J,M)*bes_cos(L,M)
            phitTMPS(J,L) = sum_top_S(J,M)*bes_sin(L,M)
         enddo
      enddo
      do J = rlwb, rupb
         phitTMP(J,:) = phitTMP(J,:) + 2.0*phitTMPC(J,:) +             &
                        2.0*phitTMPS(J,:)
      enddo
   enddo
   pott(rlwb:rupb,:) = factor * phitTMP(rlwb:rupb,:)

      if( isym == 3 ) then
         pott(1,:) = pott(2,:)
      else
         do L = 1, numphi_by_two
            pott(1,L) = pott(2,L+numphi_by_two)
            pott(1,L+numphi_by_two) = pott(2,L)
         enddo
      endif

   potp(:,numz,:) = pott
     ! done calculating top boundary potential

! if on edge of the pe grid reduce the convolution
! of G(r|r') with rho to a potential at the
! outer edge of the grid

   do L = philwb, phiupb
      do K = zlwb, zupb
         phisTMP(K,L) = sum_sid_C(K,1)
      enddo
   enddo
   do M = 2, mmax
      do L = philwb, phiupb
         do K = zlwb, zupb
            phisTMPC(K,L) = sum_sid_C(K,M)*bes_cos(L,M)
            phisTMPS(K,L) = sum_sid_S(K,M)*bes_sin(L,M)
         enddo
      enddo
      do K = zlwb, zupb
         phisTMP(K,:) = phisTMP(K,:) + 2.0*phisTMPC(K,:) +              &
                        2.0*phisTMPS(K,:)
     enddo
   enddo
   pots(zlwb:zupb,:) = factor*phisTMP(zlwb:zupb,:)
   if( (isym /= 1) ) then
      pots(1,:) = pots(2,:)
   endif
   potp(numr,:,:) = pots
  ! done calculating side boundary potential

return
end subroutine bessel
