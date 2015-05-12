subroutine binary_output(c1, c2, cc1, cc2, omsq, hm1, hm2, mass1, mass2, psi, h, &
            qfinal, initial_model_type, model_number, ra, za, phia,  &
            rb, zb, phib, rc, zc, phic, rd, zd, phid, re, ze, phie,  &
            rhm1, rhm2, rhom1, rhom2, xavg1, xavg2, separation,      &
            com, volume_factor, hem1, hem2, rhoem1, rhoem2,          &
            mass_c1, mass_c2, rho_1d, rho_c1d, rho_2e, rho_c2e,      &
            pres_d, pres_e)
  implicit none
  include 'runscf.h'
!*******************************************************************************
!
!  subroutine arguments
!

  real, dimension(maxit), intent(in) :: c1, c2, omsq, hm1, hm2, mass1, mass2, &
                                        cc1, cc2, hem1, hem2, mass_c1, mass_c2
  real, dimension(numr, numphi), intent(in) :: psi
  real, dimension(numr,numz,numphi) :: h
  integer, intent(in) :: qfinal
  integer, intent(in) :: initial_model_type
  integer, intent(in) :: ra, za, phia
  integer, intent(in) :: rb, zb, phib
  integer, intent(in) :: rc, zc, phic
  integer, intent(in) :: rd, zd, phid
  integer, intent(in) :: re, ze, phie  
  integer, intent(in) :: model_number
  real, dimension(3), intent(in) :: rhm1, rhm2
  real, intent(in) :: rhom1, rhom2
  real, intent(in) :: rhoem1, rhoem2  
  real, intent(in) :: xavg1, xavg2, separation, com
  real, intent(in) :: volume_factor
  real, intent(in) :: rho_1d, rho_c1d, rho_2e, rho_c2e, pres_d, pres_e

!
!*******************************************************************************
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

  real :: pi, grav
  common /constants/ pi, grav

!
!*******************************************************************************
!
! locall variables
!

  real, dimension(numr,numz,numphi) :: rchpot
  real, dimension(numr,numz,numphi) :: temp
  real :: gammac1, gammac2, gammae1,gammae2
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
  real :: kappac1,kappac2, kappae1,kappae2
  real :: pm1, pm2
  real :: period, omega
  real :: kepler
  real :: ret1, ret2
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
  real :: rho1i, rho2i
  integer :: louter1,  louter2
  integer :: I, J, K, L
  integer :: index
  character(len=50) :: model_template
  character(len=56) :: model_file
  integer :: phi1, phi2, phi3, phi4
  integer :: diac1, diae1, diac2, diae2, ae1, ac1, ae2, ac2
  integer, dimension(1) :: center1, center2
!
!*****************************************************************************************

  model_template = 'model_details_'

phi1 = int(numphi / 4.0) - 1
phi2 = int(numphi /  4.0) + 1
phi3 = int(3.0 * numphi / 4.0) - 1
phi4 = int(3.0 * numphi / 4.0) + 1

  gammae1 = 1.0 + 1.0/n1
  gammae2 = 1.0 + 1.0/n2
  gammac1 = 1.0 + 1.0/nc1
  gammac2 = 1.0 + 1.0/nc2

diac1=0
diae1=0
diac2=0
diae2=0
ac1=0
ae1=0
ac2=0
ae2=0

  primary = 1
  if ( mass2(qfinal) > mass1(qfinal) ) then
    primary = 2
  endif

  omega = sqrt(omsq(qfinal))

! sum up the virial pressuure
       do i = phi2, phi3
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_2e) then
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(nc2+1.0) 
               else
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(n2+1.0) 
               endif   
             enddo
          enddo
       enddo
       do i = phi4, numphi
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_1d) then
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(nc1+1.0) 
               else
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(n1+1.0)  
               endif
             enddo
          enddo
       enddo
       do i = 1, phi1
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_1d) then
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(nc1+1.0)   
               else
                 temp(k,j,i) = rhf(k)*rho(k,j,i)*h(k,j,i)/(n1+1.0)
               endif   
             enddo
          enddo
       enddo     	       
       	       
  call binary_sum(temp, ret1, ret2)
  s1 = volume_factor * ret1 
  s2 = volume_factor * ret2 
  stot = s1 + s2

! sum up the potential energy
  do K = philwb, phiupb
    do J = zlwb, zupb
      do I = rlwb, rupb
        temp(I,J,K) = rhf(I) * pot(I,J,K) * rho(I,J,K)
      enddo
    enddo
  enddo
  call binary_sum(temp, ret1, ret2)
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
  t1 = - omega * omega * volume_factor * ret1
  t2 = - omega * omega * volume_factor * ret2
  ttot = t1 + t2

  virialerr  = abs(2.0*ttot + 3.0*stot + wtot) / abs(wtot)
  virialerr1 = abs(2.0*t1   + 3.0*s1   + w1  ) / abs(w1)
  virialerr2 = abs(2.0*t2   + 3.0*s2   + w2  ) / abs(w2)

  pm1 = rhom1 * hm1(qfinal) / (nc1 + 1.0)
  pm2 = rhom2 * hm2(qfinal) / (nc2 + 1.0)
  
  kappac1 = rhom1*hm1(qfinal)/(nc1+1.0)/rhom1**(gammac1)
  kappac2 = rhom2*hm2(qfinal)/(nc2+1.0)/rhom2**(gammac2)
  kappae1 = kappac1*rho_c1d**gammac1/rho_1d**gammae1
  kappae2 = kappac2*rho_c2e**gammac2/rho_2e**gammae2

print*,"=============================================="
print*, "rho_1d", rho_1d, "rho_c1d", rho_c1d ,"FINAL" 
print*, "rho_2e", rho_2e, "rho_c2e", rho_c2e ,"FINAL" 
!print*, "nc1", nc1, "gammac1", gammac1, "gammae1", gammae1 
!print*, "nc2", nc2, "gammac2", gammac2, "gammae2", gammae2 
print*, "Qfinal", qfinal
!print*, "hm1", hm1(qfinal), "hm2", hm2(qfinal) 


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
  j1 = - 2.0 * omega * volume_factor * ret1
  j2 = - 2.0 * omega * volume_factor * ret2
  jtot = j1 + j2

! sum the internal energiies
       do i = phi2, phi3
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_2e) then
                 temp(k,j,i) = kappac2*nc2*rhf(k)*rho(k,j,i)**gammac2  
               else
                 temp(k,j,i) = kappae2*n2*rhf(k)*rho(k,j,i)**gammae2 
               endif   
             enddo
          enddo
       enddo
       do i = phi4, numphi
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_1d) then
                 temp(k,j,i) = kappac1*nc1*rhf(k)*rho(k,j,i)**gammac1   
               else
                 temp(k,j,i) = kappae1*n1*rhf(k)*rho(k,j,i)**gammae1 
               endif
             enddo
          enddo
       enddo
       do i = 1, phi1
          do j = 2, numz
             do k = 2, numr
               if (rho(k,j,i).gt.rho_1d) then
                 temp(k,j,i) = kappac1*nc1*rhf(k)*rho(k,j,i)**gammac1   
               else
                 temp(k,j,i) = kappae1*n1*rhf(k)*rho(k,j,i)**gammae1 
               endif   
             enddo
          enddo
       enddo
  call binary_sum(temp, ret1, ret2)
       e1 = volume_factor*ret1
       e2 = volume_factor*ret2
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
  vol1 = volume_factor * ret1
  vol2 = volume_factor * ret2
  reff1 = (0.75 * vol1 / pi)**(1.0/3.0)
  reff2 = (0.75 * vol2 / pi)**(1.0/3.0)

! calculate the y moment of the density distribution
  do K = philwb, phiupb
    do J = zlwb, zupb
      do I = rlwb, rupb
        temp(I,J,K) = rhf(I) * rhf(I) * sine(K) * rho(I,J,K)
      enddo
    enddo
  enddo
  call binary_sum(temp, ret1, ret2)
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
      rchpot(rlwb-1,K,L+numphi_by_two) = rho(rlwb,K,L)
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


!Find the diameter of the core and envelope in number of cells
  do i = rlwb, rupb
     if (rho(i,2,1).gt.(2*densmin)) then     
        diae1=diae1+1  
     endif
     if (rho(i,2,1).gt.(rho_1d)) then
        diac1=diac1+1
     endif
  enddo

  do i = rlwb, rupb
     if (rho(i,2,numphi/2+1).gt.(2*densmin)) then
        diae2=diae2+1
     endif
     if (rho(i,2,numphi/2+1).gt.(rho_2e)) then
        diac2=diac2+1
     endif
  enddo

!Find angular resolution (phi) of the core and the envelope 
   center1=maxloc(rho(:,2,1))
   center2=maxloc(rho(:,2,numphi/2+1)) 
!   print*, "center1 = ", center1
!   print*, "center2 = ", center2

  do i = 1, phi1
     if (rho(center1(1),2,i).gt.(2*densmin)) then
        ae1=ae1+1
     endif
     if (rho(center1(1),2,i).gt.(rho_1d)) then
        ac1=ac1+1
     endif
  enddo

  do i = phi4, numphi
     if (rho(center1(1),2,i).gt.(2*densmin)) then
        ae1=ae1+1
     endif
     if (rho(center1(1),2,i).gt.(rho_1d)) then
        ac1=ac1+1
     endif
  enddo

  do i = phi2, phi3
     if (rho(center2(1),2,i).gt.(2*densmin)) then
        ae2=ae2+1
     endif
     if (rho(center2(1),2,i).gt.(rho_2e)) then
        ac2=ac2+1
     endif
  enddo


!Write the output file
write(model_file,'(a,i6)') trim(model_template), model_number
open(unit=11, file=trim(model_file),form='formatted',status='unknown')

write(11,*)  'Model Number: ', model_number
write(11,*) 
write(11,*) 'For Star 1:'
write(11,*) 'Total Mass: ', mass1(qfinal)
write(11,*) 'Mass Core 1', mass_c1(qfinal)
write(11,*) 'Core mass ratio 1:', mass_c1(qfinal)/mass1(qfinal)
write(11,*) '(<x>,<y>): ', xavg1, yavg1
write(11,*) 'Maximum Density: ', rhom1
write(11,*) 'Core above density threshold: ', rho_1d
write(11,*) 'Polytropic Index Core: ',nc1
write(11,*) 'Polytropic Index Envelope: ',n1       
write(11,*) 'Polytropic Constant Core: ',kappac1
write(11,*) 'Polytropic Constant Envelope: ',kappae1 
write(11,*) 'Integration constant for core: ', cc1(qfinal)
write(11,*) 'Integration constant for envelope: ', c1(qfinal)
write(11,*) 'Core diameter (in number of cells):', diac1
write(11,*) 'Star diameter (in number of cells):', diae1
write(11,*) 'Core angular resolution (in number of cells):', ac1
write(11,*) 'Star angular resolution (in number of cells):', ae1
write(11,*) 'Virial Pressure: ', s1
write(11,*) 'Potential Energy: ', w1
write(11,*) 'Kinetic Energy: ', t1
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
       write(11,*) 'Roche lobe filling factor', vol1/volr1
write(11,*) 'Angular Momentum: ', j1
write(11,*) 'Internal Energy: ', e1
write(11,*) 'Total Energy: ', en1
write(11,*) 'Outer Lagrange Point: ', l2loc
write(11,*) 'Outer Edge of Roche Lobe 1: ', rochemax1
write(11,*)

write(11,*) 'For Star 2:'
write(11,*) 'Total Mass: ', mass2(qfinal)
write(11,*) 'Mass Core 2', mass_c2(qfinal)
write(11,*) 'Core mass ratio 2:', mass_c2(qfinal)/mass2(qfinal)
write(11,*) '(<x>,<y>): ', xavg2, yavg2
write(11,*) 'Maximum Density: ', rhom2
write(11,*) 'Core above density threshold: ', rho_2e
write(11,*) 'Polytropic Index Core: ',nc2
write(11,*) 'Polytropic Index Envelope: ',n2       
write(11,*) 'Polytropic Constant Core: ',kappac2
write(11,*) 'Polytropic Constant Envelope: ',kappae2 
write(11,*) 'Integration constant for core: ', cc2(qfinal)
write(11,*) 'Integration constant for envelope: ', c2(qfinal)
write(11,*) 'Core diameter (in number of cells):', diac2
write(11,*) 'Envelope diameter (in number of cells):', diae2
write(11,*) 'Core angular resolution (in number of cells):', ac2
write(11,*) 'Star angular resolution (in number of cells):', ae2
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
write(11,*) 'Roche lobe filling factor', vol2/volr2
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
write(11,*) 'Initial Model Type: ', initial_model_type
write(11,*)
write(11,*) 'rho_1d=',rho_1d,'rho_c1d=',rho_c1d
write(11,*) 'rho_2e=',rho_2e,'rho_c2e=',rho_c2e

close(11)
	 
print*, "File ", trim(model_file), " printed"

rho1i = (rho_1d+rho_c1d)/2.0
rho2i = (rho_2e+rho_c2e)/2.0

open(unit=13,file="autoread.dat")
write(13,*) nc1, " ", n1, " ", nc2, " ", n2, " ", numr, " ", numz, " ", numphi, " ",  &
     omega, " ", kappac1, " ", kappae1, " ", kappac2, " ", kappae2, " ", rho_c1d, " ",&
      rho_1d, " ",rho_c2e, " ",rho_2e, " ", pres_d, " ", pres_e, " ", xcrit, " ", com 
close(13)

end subroutine binary_output
