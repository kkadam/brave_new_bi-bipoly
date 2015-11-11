subroutine binary_scf(model_number, initial_model_type, ra, rb, rc, rd, re, rhom1, &
                      rhom2, frac, qfinal)
   implicit none
   include 'runscf.h'
!**************************************************************************************
!
!  subroutine arguments
!

integer, intent(in) :: model_number
integer, intent(in) :: initial_model_type
integer, intent(in) :: ra, rb, rc, rd, re
real, intent(in) :: rhom1
real, intent(in) :: rhom2
real, intent(in) :: frac
integer, intent(out) :: qfinal

!
!**************************************************************************************
!
!   global varaibles
!

real, dimension(numr,numz,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr) :: r, rhf, rinv, rhfinv
real, dimension(numz) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numphi) :: cosine
real, dimension(numphi) :: sine
common /trig/ cosine, sine

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numr,numz,numphi) :: pres

!
!**************************************************************************************
!
!    local variables
!

real, dimension(numr, numz, numphi) :: h, pot_it, pot_old, temp, const_map, temp_map

real, dimension(numr, numphi) :: psi

real, dimension(maxit) :: c1, c2, mass1, mass2, omsq, hm1, hm2, cc1, cc2, &
                            mass_c1, mass_c2, hem1, hem2
real :: cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2

real :: dpot, dpsi, psitmp1, psitmp2, pottmp1, pottmp2

real :: ret1, ret2

real :: global_ret1, global_ret2

real :: xavg1, xavg2, com, separation

real :: virial_error, virial_error1, virial_error2

real :: virial_error_prev

real :: volume_factor

real :: gamma

real :: time1, time2

integer :: za, phia
integer :: zb, phib
integer :: zc, phic
integer :: zd, phid
integer :: ze, phie  

real :: temp_hm1, temp_hm2

real, dimension(3) :: rhm1, rhm2, temp_rhm1, temp_rhm2

integer :: I, J, K, L, Q

  real :: pot_a, pot_b, pot_c, pot_d, pot_e
  real :: psi_a, psi_b, psi_c, psi_d, psi_e


  real :: rho_c1d, rho_1d, rho_c2e, rho_2e, pres_d, pres_e
  real :: h_c1d, h_e1d, h_c2e, h_e2e 
  real :: rhoem1, rhoem2, norm1, norm2
  integer :: rem1, zem1, phiem1, rem2, zem2, phiem2
  integer :: rm1, zm1, phim1, rm2, zm2, phim2
  integer :: rmax
  real :: gammac1, gammac2, gammae1,gammae2
  real :: kappac1,kappac2, kappae1,kappae2
  real :: rho_cc1, rho_cc2
 character*20 char1,char2,char3,char4,char5,char6,char7,char8

  integer :: diac1, diae1, diac2, diae2, ae1, ac1, ae2, ac2
  integer, dimension(1) :: center1, center2
  real :: div
  integer :: div_flag, div_it, norm_flag
  real :: x
  real, dimension(numphi) :: cos_cc
  real ::  K_part, Pi_part, W_part

!
!***********************************************************************************

   call cpu_time(time1)

! Initialize/ Renitialize arrays/ variables
   qfinal = 1
   div=0.02
   div_flag=0
   div_it=20
   norm_flag=1

   volume_factor = 2.0 * dr * dz * dphi
   rmax = numr !- 8

   h = 0.0
   psi = 0.0

   mass1 = 0.0
   mass2 = 0.0
   mass_c1=0.0
   mass_c2=0.0
   omsq = 0.0
   hm1 = 0.0
   hm2 = 0.0
   hem1 = 0.0
   hem2 = 0.0
   c1 = 0.0
   c2 = 0.0
   rho_cc1=0.0
   rho_cc2=0.0

   virial_error = 1.0

   pot_it=0
   pot_old=0
   temp=0
   const_map=0
   temp_map=0


! Initialize phi and z values of boundary points
   phia = 1
   phib = 1
   phic = numphi / 2 + 1
   phid = 1
   phie = numphi/2 + 1
   za = 2
   zb = 2
   zc = 2
   zd = 2
   ze = 2


! Find Gammas
   gammae1 = 1.0 + 1.0/n1
   gammae2 = 1.0 + 1.0/n2
   gammac1 = 1.0 + 1.0/nc1
   gammac2 = 1.0 + 1.0/nc2

! Find phi and cos(phi) arrays
   x = 0.0
   do L = 1, numphi
     phi(L) = x * dphi
     x = x + 1.0
   enddo

   do L = 1, numphi
     cos_cc(L) = cos(phi(L))
   enddo

! calculate the initial total mass
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)

mass1(1) = volume_factor * ret1
mass2(1) = volume_factor * ret2

! calculate the initial center of mass
do K = philwb, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         temp(I,J,K) = rhf(I) * rhf(I) * cosine(K) * rho(I,J,K)
      enddo
   enddo
enddo
call binary_sum(temp, ret1, ret2)

xavg1 = volume_factor * ret1 / mass1(1)
xavg2 = volume_factor * ret2 / mass2(1)
separation = xavg1 - xavg2
com = separation * mass2(1) / ( mass1(1) + mass2(1) )
com = xavg1 - com


! Open logfiles
   open(unit=13,file='iteration_log',form='formatted',status='unknown',position='append')
   write(13,*)  1, mass1(1), mass2(1), xavg1, xavg2, com, separation

   open(unit=12,file='convergence_log',form='formatted',status='unknown')

   open(unit=19,file='virial_log',form='formatted',status='unknown')



! START OF THE ITERATION CYCLE
   print*, 'maxit = ', maxit

   do Q = 2, maxit-1 
   print*, "================"
   print*, "iteration number = ", Q

   ! Solve the Poisson equation for the current density field
      if ( Q == 2 ) then
         call potential_solver(0)
         pot_old = pot
      else
         call potential_solver(1)
      endif
      pot_it = (1.0 - frac) * pot + frac * pot_old
      pot_old = pot


   ! Compute the form factor for the centrifugal potential
      do K = philwb, phiupb
         do I = rlwb-1, rupb+1
            psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
         enddo
      enddo

   ! Find pot and psi at the boundary points
      pot_a = 0.5*(pot_it(ra,za,phia) + pot_it(ra-1,za,phia))
      pot_b = 0.5*(pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib))
      pot_c = 0.5*(pot_it(rc,zc,phic) + pot_it(rc+1,zc,phic))
      pot_d = 0.5*(pot_it(rd,zd,phid) + pot_it(rd-1,zd,phid))
      pot_e = 0.5*(pot_it(re,ze,phie) + pot_it(re-1,ze,phie))
          
      psi_a = 0.5*(psi(ra,phia) + psi(ra-1,phia))          
      psi_b = 0.5*(psi(rb,phib) + psi(rb+1,phib))                    
      psi_c = 0.5*(psi(rc,phic) + psi(rc+1,phic))
      psi_d = 0.5*(psi(rd,phid) + psi(rd-1,phid))
      psi_e = 0.5*(psi(re,phie) + psi(re-1,phie))
         
!print*, pot_a, pot_b, pot_c, pot_d, pot_e
!print*, psi_a, psi_b, psi_c, psi_d, psi_e 


   ! Calculate the angular frequency and envelope c's
      omsq(q) = - (pot_a-pot_b)/(psi_a-psi_b)
   
      c1(q) = pot_b + omsq(q)*psi_b
      c2(q) = pot_c + omsq(q)*psi_c
       

   ! Calculate the core c's 
          rho_1d = 0.5*(rho(rd,zd,phid) + rho(rd-1,zd,phid))
          rho_c1d=rho_1d*muc1/mu1       
          h_e1d = c1(q) - pot_d - omsq(q)*psi_d
          h_c1d = h_e1d * (nc1+1)/(n1+1)*mu1/muc1                    
          cc1(q) = h_c1d + pot_d+ omsq(q)*psi_d 

          rho_2e = 0.5*(rho(re,ze,phie) + rho(re-1,ze,phie))
          rho_c2e=rho_2e*muc2/mu2       
          h_e2e = c2(q) - pot_e - omsq(q)*psi_e
          h_c2e = h_e2e * (nc2+1)/(n2+1)*mu2/muc2                    
          cc2(q) = h_c2e + pot_e+ omsq(q)*psi_e       


   ! Now compute the new  enthalpy field from the potential and SCF constants
          do i = 1,phi1+1
             do j = 1,numz
                do k = 1,rmax
                   if (rho(k,j,i).gt.rho_1d) then
                      h(k,j,i) = cc1(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)
                   else
                      h(k,j,i) = c1(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)  
                      
                      if(h(k,j,i).gt.hem1(q)) then
                        hem1(q) = h(k,j,i)
                        rem1 = k
                        zem1 = j
                        phiem1 = i
                      endif

                   endif
                   if(h(k,j,i).gt.hm1(q)) then
                      hm1(q) = h(k,j,i)
                      rm1 = k
                      zm1 = j
                      phim1 = i
                   endif
                enddo
             enddo
          enddo
          do i = phi2,phi3+1
             do j = 1,numz
                do k = 1,rmax
                   if (rho(k,j,i).gt.rho_2e) then
                      h(k,j,i) = cc2(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)
                   else
                      h(k,j,i) = c2(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)  
                      
                      if(h(k,j,i).gt.hem2(q)) then
                        hem2(q) = h(k,j,i)
                        rem2 = k
                        zem2 = j
                        phiem2 = i
                      endif
                      
                   endif
                   if(h(k,j,i).gt.hm2(q)) then
                      hm2(q) = h(k,j,i)
                      rm2 = k
                      zm2 = j
                      phim2 = i
                   endif
                enddo
             enddo
          enddo
          do i = phi4,numphi
             do j = 1,numz
                do k = 1,rmax
                   if (rho(k,j,i).gt.rho_1d) then
                      h(k,j,i) = cc1(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)
                   else
                      h(k,j,i) = c1(q) - pot_it(k,j,i) - omsq(q)*psi(k,i)  
                      
                      if(h(k,j,i).gt.hem1(q)) then
                        hem1(q) = h(k,j,i)
                        rem1 = k
                        zem1 = j
                        phiem1 = i
                      endif
                      
                   endif
                   if(h(k,j,i).gt.hm1(q)) then
                      hm1(q) = h(k,j,i)
                      rm1 = k
                      zm1 = j
                      phim1 = i
                   endif
                enddo
             enddo
          enddo 

          do i = 1,numphi
             do j = 1,numz
                do k = 1,rmax
                    if (h(k,j,i).lt.0) then
                       h(k,j,i)=0.0
                    endif
                enddo
             enddo
          enddo

   
   ! Calculate normalization constants
          norm1 = mu1/muc1*(h_c1d/hm1(q))**nc1        
          norm2 = mu2/muc2*(h_c2e/hm2(q))**nc2

   ! Calculate the new density field from the enthalpy
      rho_cc1=0.0
      rho_cc2=0.0

          do i = 1,phi1+1
             do j = 1,numz
                do k = 1,numr
                   if(h(k,j,i).gt.0.0) then
                      if (rho(k,j,i).gt.rho_1d) then
                         rho(k,j,i) = rhom1*(h(k,j,i)/hm1(q))**nc1
                      else
                         !rho(k,j,i) = rho_c1d*mu1/muc1*(h(k,j,i)/h_e1d)**n1
                         rho(k,j,i) = rhom1*norm1*(h(k,j,i)/h_e1d)**n1
                      endif
                   else   
                      rho(k,j,i) = 0.0
                   endif

                   if(rho(k,j,i).gt.rho_cc1) then
                      rho_cc1 = rho(k,j,i)
                   endif

                enddo
             enddo
          enddo
          do i = phi2,phi3+1
             do j = 1,numz
                do k = 1,numr
                   if(h(k,j,i).gt.0.0) then
                      if (rho(k,j,i).gt.rho_2e) then
                         rho(k,j,i) = rhom2*(h(k,j,i)/hm2(q))**nc2
                      else
                         !rho(k,j,i) = rho_c2e*mu2/muc2*(h(k,j,i)/h_e2e)**n2
                         rho(k,j,i) = rhom2*norm2*(h(k,j,i)/h_e2e)**n2
                      endif
                   else   
                      rho(k,j,i) = 0.0
                   endif

                   if(rho(k,j,i).gt.rho_cc2) then
                      rho_cc2 = rho(k,j,i)
                   endif

                enddo
             enddo
          enddo
          do i = phi4,numphi
             do j = 1,numz
                do k = 1,numr
                   if(h(k,j,i).gt.0.0) then
                      if (rho(k,j,i).gt.rho_1d) then
                         rho(k,j,i) = rhom1*(h(k,j,i)/hm1(q))**nc1
                      else
                         !rho(k,j,i) = rho_c1d*mu1/muc1*(h(k,j,i)/h_e1d)**n1
                         rho(k,j,i) = rhom1*norm1*(h(k,j,i)/h_e1d)**n1
                      endif
                   else   
                      rho(k,j,i) = 0.0
                   endif

                   if(rho(k,j,i).gt.rho_cc1) then
                      rho_cc1 = rho(k,j,i)
                   endif

                enddo
             enddo
          enddo

   ! Conditionally normalize wrt the central densities of the stars
      if ( norm_flag == 1 ) then
          print*, "normalizing "
          do i = 1,phi1+1
             do j = 2,numz
                do k = 2,numr
                   rho(k,j,i) = rho(k,j,i)/rho_cc1*rhom1
                enddo
             enddo
          enddo

          do i = phi2,phi3+1
             do j = 2,numz
                do k = 2,numr
                   rho(k,j,i) = rho(k,j,i)/rho_cc2*rhom2
                enddo
             enddo
          enddo


          do i = phi4,numphi
             do j = 2,numz
                do k = 2,numr
                   rho(k,j,i) = rho(k,j,i)/rho_cc1*rhom1
                enddo
             enddo
          enddo
      endif

!print*, rm1,zm1,phim1
!print*, rm2,zm2,phim2


   ! Zero out the density field between the axis and the inner boundary points
   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         do J = rlwb-1, rupb
            if( ( rhf(J) * cos_cc(L) .ge. 0.0  ) .and. &
                ( rhf(J) * cos_cc(L) .lt. rhf(rb) * 1.0 ) ) then
                 rho(J,K,L) = 0.0
            endif
          enddo
       enddo
    enddo

   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         do J = rlwb-1, rupb
            if( ( rhf(J) * cos_cc(L) .le. 0.0  ) .and. &
                ( rhf(J) * cos_cc(L) .gt. rhf(rc) * -1.0 ) ) then
                 rho(J,K,L) = 0.0
            endif
          enddo
       enddo
    enddo


   ! zero out the density field on RHS of point A
   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         do J = rlwb-1, rupb
            if( ( rhf(J) * cos_cc(L) .gt. rhf(ra)) ) then
                 rho(J,K,L) = 0.0
            endif
          enddo
       enddo
    enddo


   ! impose the equatorial boundary condition
      do K = philwb, phiupb
         do I = rlwb, rupb
            rho(I,zlwb-1,K) = rho(I,zlwb,K)
         enddo
      enddo

   ! impose the axial boundary condition
      rho(rlwb-1,:,:) = cshift(rho(rlwb,:,:),dim=2,shift=numphi/2)


   ! has the solution converged?
   cnvgom = abs( (omsq(Q) - omsq(Q-1)) / omsq(Q) )
   cnvgc1 = abs( (c1(Q) - c1(Q-1)) / c1(Q) )
   cnvgc2 = abs( (c2(Q) - c2(Q-1)) / c2(Q) )
   cnvgh1 = abs( (hm1(Q) - hm1(Q-1)) / hm1(Q) )
   cnvgh2 = abs( (hm2(Q) - hm2(Q-1)) / hm2(Q) )

   virial_error_prev = virial_error

   call compute_virial_field(psi, rho_1d, rho_2e, h, sqrt(omsq(Q)), volume_factor, &
                             virial_error1, virial_error2, virial_error, K_part, Pi_part, W_part)

! Calculating stuff for printing >>


   ! calculate the total mass for each star  and core mass for each star
   do K = philwb, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)

   mass1(Q) = volume_factor * ret1 
   mass2(Q) = volume_factor * ret2

          temp = 0.0          
          do i = 1, numphi
             do j = 2, numz
                do k = 2, numr
                   if (rho(k,j,i).gt.rho_1d) then 
                     temp(k,j,i) = rhf(k)*rho(k,j,i)
                   endif
                enddo
             enddo
          enddo          
          call binary_sum(temp, ret1, ret2)
          mass_c1(q) = volume_factor*ret1
          
          temp=0.0
          do i = 1, numphi
             do j = 2, numz
                do k = 2, numr
                   if (rho(k,j,i).gt.rho_2e) then 
                     temp(k,j,i) = rhf(k)*rho(k,j,i)
                   endif
                enddo
             enddo
          enddo          
          call binary_sum(temp, ret1, ret2)
          mass_c2(q) = volume_factor*ret2  

!print*, "mass1  = ", mass1(Q)
!print*, "massc1 = ", mass_c1(Q) 
!print*, "mass2  = ", mass2(Q)
!print*, "massc2 = ", mass_c2(Q)

   ! calculate the center of mass for each star
   do K = philwb, phiupb
      do J  = zlwb, zupb
         do I = rlwb, rupb
            temp(I,J,K) = rhf(I) * rhf(I) * cosine(K) * rho(I,J,K)
         enddo
      enddo
   enddo
   call binary_sum(temp, ret1, ret2)

   xavg1 = volume_factor * ret1 / mass1(Q)
   xavg2 = volume_factor * ret2 / mass2(Q)
   separation = xavg1 - xavg2
   com = separation * mass2(Q) / (mass1(Q) + mass2(Q) )
   com = xavg1 - com
 
!print*, "omsq = ", cnvgom
!print*, "c1 = ",  cnvgc1,  "c2 = ",  cnvgc2
!print*,  "h1 = ",  cnvgh1,  "h2 = ",  cnvgh2
!print*, "hm1(q) = ",hm1(q),"hm2(q) = ",hm2(q)
!print*, "rmom1 = ", rhom1, "rhom2 = ", rhom2
!print*, "",,"",


!Find the diameter of the core and envelope in number of cells
  diac1=0
  diac2=0
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

  ac1=0
  ac2=0
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

write (char1, "(F10.7)") cnvgom
write (char2, "(F10.7)") cnvgc1
write (char3, "(F10.7)") cnvgc2
write (char4, "(F10.7)") cnvgh1
write (char5, "(F10.7)") cnvgh2
write (char6, "(F10.7)") rho_cc1
write (char7, "(F10.7)") rho_cc2
write (char8, "(F10.7)") virial_error

       print*, "mass ratio = ", "m1/m2=",mass1(Q)/mass2(Q),"or m2/m1=",mass2(Q)/mass1(Q)
       print*, "core mass ratio 1 = ", mass_c1(Q)/mass1(Q)
       print*, "core mass ratio 2 = ", mass_c2(Q)/mass2(Q)
       print*, "core1 resolution: ", diac1, ac1 
       print*, "core2 resolution: ", diac2, ac2

      write(13,*) rho_cc1, rho_cc2,hm1(Q),hm2(Q)

      write(12,*) trim(char1),trim(char2),trim(char3),trim(char4),trim(char5),trim(char8)
      write(19,*) K_part, Pi_part, W_part 
! Finished printing stuff ^^

! Convergance test
   if ( cnvgom < eps .and. cnvgc1 < eps .and. cnvgc2 < eps .and. &
        cnvgh1 < eps .and. cnvgh2 < eps ) then!.and. virial_error < eps ) then
      exit
   elseif ( Q > div_it .and. (cnvgom > div .or. cnvgc1 > div .or. cnvgc2 > div .or. &
        cnvgh1  > div .or. cnvgh2 > div)) then
      print*,"====================================="
      print*,"####   Solution diverged. :(    ####"
      print*,"====================================="
      div_flag=1
      exit
   elseif (mass1(Q)/mass2(Q) > div_it .or. mass2(Q)/mass1(Q) > div_it) then
      print*,"====================================================="
      print*,"####    Divergence by extreme mass ratio. X(    ####"
      print*,"====================================================="
      div_flag=1
      exit   
   endif 

!   if ( virial_error > virial_error_prev .and. Q > 10  ) then
!      exit
!   endif


enddo                   
! END OF THE ITERATION CYCLE

print*,"SEMIFINAL Q =", Q
if ( q >= maxit ) then
   qfinal = maxit
else
   qfinal = q + 1
endif

print*,"FINAL Q =", Q

if (q==maxit .and. ( cnvgom > eps .or. cnvgc1 > eps .or. cnvgc2 > eps .or. &
        cnvgh1 > eps .or. cnvgh2 > eps) ) then
   print*,"SATISFIED"
   div_flag=1
endif

  rhm1(1) = rhf(rm1)
  rhm1(2) = zhf(zm1)
  rhm1(3) = phi(phim1)

!  rhm2(1) = rhf(rm1)      !wtheck?
  rhm2(1) = rhf(rm2) 
  rhm2(2) = zhf(zm2)
  rhm2(3) = phi(phim2)

! syncchronize the potential with the converged density distribution
call potential_solver(1)

! and synchronize the centrifugal potential
do K = philwb, phiupb
   do I = rlwb-1, rupb+1
      psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
   enddo
enddo

! now calculate the final angular frequency
          pot_a = 0.5*(pot_it(ra,za,phia) + pot_it(ra-1,za,phia))
          pot_b = 0.5*(pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib))
          pot_c = 0.5*(pot_it(rc,zc,phic) + pot_it(rc+1,zc,phic))
          pot_d = pot_it(rd,zd,phid)
          pot_e = pot_it(re,ze,phie)
          
          psi_a = 0.5*(psi(ra,phia) + psi(ra-1,phia))          
          psi_b = 0.5*(psi(rb,phib) + psi(rb+1,phib))                    
          psi_c = 0.5*(psi(rc,phic) + psi(rc+1,phic))
          psi_d = psi(rd,phid)
          psi_e = psi(re,phie)
          
          omsq(qfinal) = - (pot_a-pot_b)/(psi_a-psi_b) 

! and calculate the two integration constants
          c1(qfinal) = pot_b + omsq(q)*psi_b
          c2(qfinal) = pot_c + omsq(q)*psi_c
       
          rho_1d = rho(rd,zd,phid)
          rho_c1d= rho_1d*muc1/mu1       
          h_e1d = c1(q) - pot_d - omsq(q)*psi_d
          h_c1d = h_e1d * (nc1+1)/(n1+1)*mu1/muc1                    
          cc1(qfinal) = h_c1d + pot_d+ omsq(q)*psi_d 
        
          rho_2e = rho(re,ze,phie)
          rho_c2e = rho_2e*muc2/mu2       
          h_e2e = c2(q) - pot_e - omsq(q)*psi_e
          h_c2e = h_e2e * (nc2+1)/(n2+1)*mu2/muc2                    
          cc2(qfinal) = h_c2e + pot_e+ omsq(q)*psi_e


  hm1(qfinal) = hm1(qfinal-1)
  hm2(qfinal) = hm2(qfinal-1)
  hem1(qfinal) = hem1(qfinal-1)
  hem2(qfinal) = hem2(qfinal-1)  
  mass1(qfinal) = mass1(qfinal-1)
  mass2(qfinal) = mass2(qfinal-1)
  mass_c1(qfinal) = mass_c1(qfinal-1)
  mass_c2(qfinal) = mass_c2(qfinal-1)
  rhoem1 = rho(rem1,zem1,phiem1)
  rhoem2 = rho(rem2,zem2,phiem2)  


! Calculate pressure
!  print*,"hm1",hm1(qfinal)
   kappac1 = rhom1*hm1(qfinal)/(nc1+1.0)/rhom1**(gammac1)
   kappac2 = rhom2*hm2(qfinal)/(nc2+1.0)/rhom2**(gammac2)
   kappae1 = kappac1*rho_c1d**gammac1/rho_1d**gammae1
   kappae2 = kappac2*rho_c2e**gammac2/rho_2e**gammae2

!  print*, "kappas used for pressure file:"
!  print*, "kappac1=",kappac1,"kappae1=",kappae1
!  print*, "kappac2=",kappac2,"kappae2=",kappae2

   call compute_pressure(rho,pres,kappac1,kappae1,kappac2,kappae2,rho_1d,rho_c1d,rho_2e,rho_c2e)
          pres_d = pres(rd,zd,phid)
          pres_e = pres(re,ze,phie)
 
 
!   write(13,*) 'Model: ', model_number, ' done in time: ', time2 - time1

   call binary_output(c1, c2, cc1, cc2, omsq, hm1, hm2, mass1, mass2, psi, h, &
            qfinal, initial_model_type, model_number, ra, za, phia,          &
            rb, zb, phib, rc, zc, phic, rd, zd, phid, re, ze, phie,          &
            rhm1, rhm2, rho_cc1, rho_cc2, xavg1, xavg2, separation,              &
            com, volume_factor, hem1, hem2, rhoem1, rhoem2,                  &
            mass_c1, mass_c2, rho_1d, rho_c1d, rho_2e, rho_c2e,              &
            pres_d, pres_e, rem1, rem2, div_flag)


!  call ancient_output(c1, c2, omsq, hm1, hm2, mass1, mass2, psi, h, qfinal,  &
!                   initial_model_type, model_number, ra, za, phia, rb, zb,   &
!                   phib, rc, zc, phic, rhm1, rhm2, 1.5, rhom1, rhom2, xavg1, &
!                   xavg2, separation, com, volume_factor)



!call cpu_time(time1)


! Print various maps
    do i = 1, numphi
       do j = 1, numz
          do k = 1, numr
             const_map(k,j,i)=h(k,j,i)+pot(k,j,i)+omsq(qfinal)*psi(k,j)
          enddo
       enddo
    enddo

    do i = 1, numphi
       do j = 1, numz
          do k = 1, numr
             temp_map(k,j,i)=pot(k,j,i)+omsq(qfinal)*psi(k,j)
          enddo
       enddo
    enddo

   call output('density.bin','star',rho)
   call output('pressure.bin','pres',pres)
   call output('const_map.bin','const',const_map)
   call output('enthalpy.bin','enth',h)
   call output('pot.bin','pot',pot)
   call output('pot_eff.bin','eff',temp_map)

   open(unit=10,file='psi')
      do j=1,numz
         do i=1,numr
            write(10,*) i,j,psi(i,j)
         enddo
         write(10,*)
      enddo
   close(10)
   print*, "File psi printed"

!call cpu_time(time2)

!   write(13,*) 'Model: ', model_number, ' disk I/O done in time: ', time2 - time1
   close(19)
   close(13)
   close(12)

end subroutine binary_scf
