subroutine binary_initialize(rhom1, rhom2, ra, rb, rc, phia, phib, phic, &
                             za, zb, zc, initial_model_type)
implicit none
include 'runscf.h'
!****************************************************************************************
!
!  subroutine arguments
!

real, intent(in) :: rhom1
real, intent(in) :: rhom2
integer, intent(in) :: ra
integer, intent(in) :: rb
integer, intent(in) :: rc
integer, intent(in) :: phia
integer, intent(in) :: phib
integer, intent(in) :: phic
integer, intent(in) :: za
integer, intent(in) :: zb
integer, intent(in) :: zc
integer, intent(in) :: initial_model_type

!
!  initial_model_type = 1 ==> spherical Gausssian
!                     = 2 ==> uniform density sphere
!                     = 3 ==> n=1 polytrope
!                     = 4 ==> double gaussian
!                     = 5 ==> double gaussian with density jump
!                     = 6 ==> use existing density, file nnamed density.init
!
!****************************************************************************************
!
!  global variables
!

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr) :: rhf, r, rhfinv, rinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numphi) :: cosine, sine
common /trig/ cosine, sine


!
!****************************************************************************************
!
!  local variables
!

integer :: rindex_center_1, rindex_center_2

real :: sigsq1, sigsq2, sigsqc1, sigsqc2

real :: radius_star1, radius_star2

real :: k1, k2

real :: dsq, d

integer :: I, J, K

real :: th1, th2, thc1, thc2,avth1, avth2
!
!****************************************************************************************

select case(initial_model_type)

   case(1)
   ! gaussian density field
   rindex_center_1 = (ra - rb) / 2 + rb
   rindex_center_2 = (ra - rc) / 2 + rc

   dsq = (rhf_g(rindex_center_1) - rhf_g(rb)) * (rhf_g(rindex_center_1) - rhf_g(rb))
   sigsq1 = -dsq / log(1.0e-6)

   dsq = (rhf_g(rindex_center_2) - rhf_g(rc)) * (rhf_g(rindex_center_2) - rhf_g(rc))
   sigsq2 = - dsq / log(1.0e-6)

   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) *  &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) +  &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia))*       &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +      &
                  zhf(J)*zhf(J)
            rho(I,J,K) = rhom1 * exp(-dsq / sigsq1)
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) *  &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) +  &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) *      &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) +      &
                  zhf(J)*zhf(J)
            rho(I,J,K) = rhom2 * exp(-dsq / sigsq2)
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +     &
                  zhf(J)*zhf(J)
            rho(I,J,K) = rhom1 * exp(-dsq / sigsq1)
         enddo
      enddo
   enddo

   case(2)
   ! uniform density sphere
   rindex_center_1 = (ra - rb) / 2 + rb
   rindex_center_2 = (ra - rc) / 2 + rc
   radius_star1 = (rhf_g(rindex_center_1) - rhf_g(rb))*(rhf_g(rindex_center_1) - rhf_g(rb))
   radius_star2 = (rhf_g(rindex_center_2) - rhf_g(rc))*(rhf_g(rindex_center_2) - rhf_g(rc))
   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +     &
                  zhf(J)*zhf(J)
            if ( dsq <= radius_star1 ) then
               rho(I,J,K) = rhom1
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) +     &
                  zhf(J)*zhf(J)
            if ( dsq <= radius_star2 ) then
               rho(I,J,K) = rhom2
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +     &
                  zhf(J)*zhf(J)
            if ( dsq <- radius_star1 ) then
               rho(I,J,K) = rhom1
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   case(3)
   ! n=1 polytrope with rho = sinc(r) = sin(r) / r
   rindex_center_1 = (ra - rb) / 2 + rb
   rindex_center_2 = (ra - rc) / 2 + rc

   radius_star1 = (rhf_g(rindex_center_1) - rhf_g(rb))*(rhf_g(rindex_center_1) - rhf_g(rb))
   radius_star1 = sqrt(radius_star1)
   k1 = pi / radius_star1
   radius_star2 = (rhf_g(rindex_center_2) - rhf_g(rc))*(rhf_g(rindex_center_2) - rhf_g(rc))
   radius_star2 = sqrt(radius_star2)
   k2 = pi / radius_star2
   do K = philwb, phi1
      do  J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +     &
                  zhf(J)*zhf(J)
            d = sqrt(dsq)
            if ( d <= radius_star1 ) then
               if ( d > 0.0 ) then
                  rho(I,J,K) = rhom1 * sin(k1 * d) / k1 / d
               else
                  rho(I,J,K) = rhom1
               endif
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) +     &
                  zhf(J)*zhf(J)
            d = sqrt(dsq)
            if ( d <= radius_star2) then
               if ( d > 0.0 ) then
                  rho(I,J,K) = rhom2 * sin(k2 * d) / k2 / d
               else
                  rho(I,J,K) = rhom2
               endif
            else
                rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +    &
                  zhf(J)*zhf(J)
            d = sqrt(dsq)
            if ( d <= radius_star1 ) then
               if ( d > 0.0 ) then
                  rho(I,J,K) = rhom1 * sin(k1 * d) / k1 / d
               else
                  rho(I,J,K) = rhom1
               endif
            else
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   case(4)
   ! Double gaussian model
   rindex_center_1 = (ra - rb) / 2 + rb
   rindex_center_2 = (ra - rc) / 2 + rc

   dsq = (rhf_g(rindex_center_1) - rhf_g(rb)) * (rhf_g(rindex_center_1) - rhf_g(rb))
   sigsq1 = -dsq / log(1.0e-2)
   sigsqc1 = -dsq / log(1.0e-11)

   dsq = (rhf_g(rindex_center_2) - rhf_g(rc)) * (rhf_g(rindex_center_2) - rhf_g(rc))
   sigsq2 = - dsq / log(1.0e-2)
   sigsqc2 = -dsq / log(1.0e-11)

   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) *  &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) +  &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia))*       &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +      &
                  zhf(J)*zhf(J)
            rho(I,J,K) = max(rhom1/4 * exp(-dsq / sigsq1), rhom1 * exp(-dsq / sigsqc1)) 
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) *  &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) +  &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) *      &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) +      &
                  zhf(J)*zhf(J)
            rho(I,J,K) = max(rhom2/4 * exp(-dsq / sigsq2),rhom2 * exp(-dsq / sigsqc2)) 
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +     &
                  zhf(J)*zhf(J)
            rho(I,J,K) = max(rhom1/4 * exp(-dsq / sigsq1), rhom1 * exp(-dsq / sigsqc1))
         enddo
      enddo
   enddo

   do K = philwb, phiupb
      do J = zlwb-1, zupb+1
         do I = rlwb-1, rupb+1
             if (rho(I,J,K).lt.1e-3) then
                 rho(I,J,K) = 0.0
             endif
         enddo
      enddo
   enddo

   case(5)
   ! Double gaussian model with density gap
   rindex_center_1 = (ra - rb) / 2 + rb
   rindex_center_2 = (ra - rc) / 2 + rc

   dsq = (rhf_g(rindex_center_1) - rhf_g(rb)) * (rhf_g(rindex_center_1) - rhf_g(rb))
   sigsq1 = -dsq / log(1.0e-2)
   sigsqc1 = -dsq / log(1.0e-11)

   dsq = (rhf_g(rindex_center_2) - rhf_g(rc)) * (rhf_g(rindex_center_2) - rhf_g(rc))
   sigsq2 = - dsq / log(1.0e-2)
   sigsqc2 = -dsq / log(1.0e-11)

   thc1=rhom1/4*1.6
   th1=rhom1/4.5
   avth1=(thc1+th1)/2

   thc2=rhom2/4*1.6
   th2=rhom2/4.5
   avth2=(thc2+th2)/2

   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) *  &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) +  &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia))*       &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +      &
                  zhf(J)*zhf(J)
            rho(I,J,K) = max(rhom1/4 * exp(-dsq / sigsq1), rhom1 * exp(-dsq / sigsqc1))
            if ( (rho(I,J,K).lt.thc1).and.(rho(I,J,K).gt.avth1) ) then
                rho(I,J,K) = thc1
            elseif ( (rho(I,J,K).lt.avth1).and.(rho(I,J,K).gt.th1) ) then
                rho(I,J,K) = th1
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) *  &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_2)*cosine(phic)) +  &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) *      &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_2)*sine(phic)) +      &
                  zhf(J)*zhf(J)
            rho(I,J,K) = max(rhom2/4 * exp(-dsq / sigsq2),rhom2 * exp(-dsq / sigsqc2))
            if ( (rho(I,J,K).lt.thc2).and.(rho(I,J,K).gt.avth2) ) then
                rho(I,J,K) = thc2
            elseif ( (rho(I,J,K).lt.avth2).and.(rho(I,J,K).gt.th2) ) then
                rho(I,J,K) = th2
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            dsq = (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) * &
                  (rhf(I)*cosine(K) - rhf_g(rindex_center_1)*cosine(phia)) + &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) *     &
                  (rhf(I)*sine(K) - rhf_g(rindex_center_1)*sine(phia)) +     &
                  zhf(J)*zhf(J)
            rho(I,J,K) = max(rhom1/4 * exp(-dsq / sigsq1), rhom1 * exp(-dsq / sigsqc1))
            if ( (rho(I,J,K).lt.thc1).and.(rho(I,J,K).gt.avth1) ) then
                rho(I,J,K) = thc1
            elseif ( (rho(I,J,K).lt.avth1).and.(rho(I,J,K).gt.th1) ) then
                rho(I,J,K) = th1
            endif
         enddo
      enddo
   enddo

   do K = philwb, phiupb
      do J = zlwb-1, zupb+1
         do I = rlwb-1, rupb+1
             if (rho(I,J,K).lt.1e-3) then
                 rho(I,J,K) = 0.0
             endif
         enddo
      enddo
   enddo


   case(6)
   ! read in a density model from the file density.init
   open(unit=50,file='density.init',form='unformatted',convert='BIG_ENDIAN',status='unknown')
   read(50) rho
   close(50)

end select

! impose equatorial boundary condition
   rho(:,zlwb-1,:) = rho(:,zlwb,:)

! impose the axial boundary condition
   rho(rlwb-1,:,:) = cshift(rho(rlwb,:,:),dim=2,shift=numphi/2)


!   call output('rhoinit.bin','rhoinit',rho)


end subroutine binary_initialize
