subroutine binary_scf(model_number, initial_model_type, ra, rb, rc, rd, re, &
                      rhom1, rhom2, frac, eps, qfinal)
  implicit none
  include 'runscf.h'
!**************************************************************************************************
!
!  subroutine arguments
!

  integer, intent(in) :: model_number
  integer, intent(in) :: initial_model_type
  integer, intent(in) :: ra, rb, rc, rd, re
  real, intent(in) :: rhom1, rhom2
  real, intent(in) :: frac, eps
  integer, intent(out) :: qfinal

!
!**************************************************************************************************
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

!
!**************************************************************************************************
!
!    local variables
!

  real, dimension(numr, numz, numphi) :: h, pot_it, pot_old, temp
  real, dimension(numr, numphi) :: psi
  real, dimension(maxit) :: c1, c2, mass1, mass2, omsq, hm1, hm2, cc1, cc2, &
                            mass_c1, mass_c2, hem1, hem2
  real :: cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2
  real :: dpot, dpsi, psitmp1, psitmp2, pottmp1, pottmp2
  real :: ret1, ret2
  real :: xavg1, xavg2, com, separation
  real :: virial_error, virial_error1, virial_error2
  real :: virial_error_prev
  real :: volume_factor
  real :: time1, time2
  integer :: za, phia
  integer :: zb, phib
  integer :: zc, phic
  integer :: zd, phid
  integer :: ze, phie  
  real, dimension(3) :: rhm1, rhm2

  real :: pot_a, pot_b, pot_c, pot_d, pot_e
  real :: psi_a, psi_b, psi_c, psi_d, psi_e
  integer :: I, J, K, L, Q

  real :: rho_c1d, rho_1d, rho_c2e, rho_2e
  real :: h_c1d, h_e1d, h_c2e, h_e2e 
  real :: rhoem1, rhoem2, norm1, norm2
  real :: rem1, zem1, phiem1, rem2, zem2, phiem2
  real :: rm1, zm1, phim1, rm2, zm2, phim2
  integer :: rmax

!
!**************************************************************************************************

  call cpu_time(time1)

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
  
  volume_factor = 2.0 * dr * dz * dphi

  rmax = numr - 8

!  Re-Initialize arrays/ variables  
  h = 0.0
  psi = 0.0

  qfinal = 1  
  
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

  virial_error = 1.0

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

  open(unit=13,file='iteration_log',form='formatted',status='unknown',position='append')
  write(13,*) mass1(1), mass2(1), xavg1, xavg2, com, separation

  
! START OF THE ITERATION CYCLE  
  do Q = 2, maxit-1                                   

   ! solve the Poisson equation for the current density field
    if ( Q == 2 ) then
      call potential_solver(0)
      do L = 1, numphi
        do K = 1, numz
          do J = 1, numr
            pot_old(J,K,L) = pot(J,K,L)
          enddo
        enddo
      enddo
    else
      call potential_solver(1)
    endif
   
    do L = 1, numphi
      do K = 1, numz
        do J = 1, numr
           pot_it(J,K,L) = (1.0 - frac) * pot(J,K,L) + frac * pot_old(J,K,L) !!wth?
           pot_old(J,K,L) = pot(J,K,L)
        enddo
      enddo
    enddo

   ! compute the form factor for the centrifugal potential
    do K = philwb, phiupb
      do I = rlwb-1, rupb+1
        psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
      enddo
    enddo

          pot_a = 0.5*(pot_it(ra,za,phia) + pot_it(ra-1,za,phia))
          pot_b = 0.5*(pot_it(rb,zb,phib) + pot_it(rb+1,zb,phib))
          pot_c = 0.5*(pot_it(rc,zc,phic) + pot_it(rc+1,zc,phic))
          pot_d = 0.5*(pot_it(rd,zd,phid) + pot_it(rd+1,zd,phid))
          pot_e = 0.5*(pot_it(re,ze,phie) + pot_it(re+1,ze,phie))
          
          psi_a = 0.5*(psi(ra,phia) + psi(ra-1,phia))          
          psi_b = 0.5*(psi(rb,phib) + psi(rb+1,phib))                    
          psi_c = 0.5*(psi(rc,phic) + psi(rc+1,phic))
          psi_d = 0.5*(psi(rd,phid) + psi(rd+1,phid))
          psi_e = 0.5*(psi(re,phie) + psi(re+1,phie))
          

          omsq(q) = - (pot_a-pot_b)/(psi_a-psi_b)
          
          c1(q) = pot_b + omsq(q)*psi_b
          c2(q) = pot_c + omsq(q)*psi_c
       
          rho_1d = 0.5*(rho(rd,zd,phid) + rho(rd+1,zd,phid))
          rho_c1d=rho_1d*muc1/mu1       
          h_e1d = c1(q) - pot_d - omsq(q)*psi_d
          h_c1d = h_e1d * (nc1+1)/(n1+1)*mu1/muc1                    
          cc1(q) = h_c1d + pot_d+ omsq(q)*psi_d 
          
          rho_2e = 0.5*(rho(re,ze,phie) + rho(re+1,ze,phie))
          rho_c2e=rho_2e*muc2/mu2       
          h_e2e = c2(q) - pot_e - omsq(q)*psi_e
          h_c2e = h_e2e * (nc2+1)/(n2+1)*mu2/muc2                    
          cc2(q) = h_c2e + pot_e+ omsq(q)*psi_e           


!  Calculate new enthalpy field and max enthalpies and their positions

          do i = 1,phi1+1
             do j = 2,numz
                do k = 2,rmax
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
             do j = 2,numz
                do k = 2,rmax
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
             do j = 2,numz
                do k = 2,rmax
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
   
          
          norm1 = mu1/muc1*(h_c1d/hm1(q))**nc1        
          norm2 = mu2/muc2*(h_c2e/hm2(q))**nc2
          

   ! calculate the new density field from the enthalpy
          do i = 1,phi1+1
             do j = 2,numz
                do k = 2,numr
                   if(h(k,j,i).gt.0.0) then
                      if (rho(k,j,i).gt.rho_1d) then	   
                         rho(k,j,i) = rhom1*(h(k,j,i)/hm1(q))**nc1
                      else
                      	 rho(k,j,i) = rhom1*norm1*(h(k,j,i)/h_e1d)**n1
                         !rho(k,j,i) = rhom1*norm1*(h(k,j,i)/hm1(q))**n1
                      endif
                   else   
                      rho(k,j,i) = 0.0	   
                   endif
                enddo
             enddo
          enddo
          do i = phi2,phi3+1
             do j = 2,numz
                do k = 2,numr
                   if(h(k,j,i).gt.0.0) then
                      if (rho(k,j,i).gt.rho_2e) then	   
                         rho(k,j,i) = rhom2*(h(k,j,i)/hm2(q))**nc2
                      else
                      	 rho(k,j,i) = rhom2*norm2*(h(k,j,i)/h_e2e)**n2
                         !rho(k,j,i) = rhom2*norm2*(h(k,j,i)/hm2(q))**n2
                      endif
                   else   
                      rho(k,j,i) = 0.0	   
                   endif
                enddo
             enddo
          enddo
          do i = phi4,numphi
             do j = 2,numz
                do k = 2,numr
                   if(h(k,j,i).gt.0.0) then
                      if (rho(k,j,i).gt.rho_1d) then	   
                         rho(k,j,i) = rhom1*(h(k,j,i)/hm1(q))**nc1
                      else
                      	 rho(k,j,i) = rhom1*norm1*(h(k,j,i)/h_e1d)**n1
                         !rho(k,j,i) = rhom1*norm1*(h(k,j,i)/hm1(q))**n1
                      endif
                   else   
                      rho(k,j,i) = 0.0	   
                   endif
                enddo
             enddo
          enddo

   ! zero out the density field between the axis and the inner boundary points
   do K = philwb, phi1
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf(rc) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   ! impose the equatorial boundary condition
   do K = philwb, phiupb
      do J = rlwb, rupb
         rho(J,zlwb-1,K) = rho(J,zlwb,K)
      enddo
   enddo

   ! impose the axial boundary condition
   do L = 1, numphi_by_two
      do K = zlwb, zupb
         rho(rlwb-1,K,L)               = rho(rlwb,K,L+numphi_by_two)
         rho(rlwb-1,K,L+numphi_by_two) = rho(rlwb,K,L)
      enddo
   enddo

   ! calculate the total mass and core mass for each star
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
 
   ! has the solution converged?
   cnvgom = abs( (omsq(Q) - omsq(Q-1)) / omsq(Q) )
   cnvgc1 = abs( (c1(Q) - c1(Q-1)) / c1(Q) )
   cnvgc2 = abs( (c2(Q) - c2(Q-1)) / c2(Q) )
   cnvgh1 = abs( (hm1(Q) - hm1(Q-1)) / hm1(Q) )
   cnvgh2 = abs( (hm2(Q) - hm2(Q-1)) / hm2(Q) )

   virial_error_prev = virial_error
   call compute_virial_error(psi, h, sqrt(omsq(Q)), rho_1d, rho_2e, &
                     volume_factor, virial_error1, virial_error2, virial_error)


   write(13,*) Q, mass1(Q), mass2(Q), xavg1, xavg2, com, omsq(Q), c1(Q), c2(Q), &
               hm1(Q), hm2(Q), cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2, virial_error1, &
               virial_error2, virial_error

   if ( cnvgom < eps .and. cnvgc1 < eps .and. cnvgc2 < eps .and. &
        cnvgh1 < eps .and. cnvgh2 < eps ) then
      exit
   endif

!   if ( virial_error > virial_error_prev .and. Q > 10  ) then
!      exit
!   endif

enddo                                               
! END OF THE ITERATION CYCLE

	
if ( q >= maxit ) then
   qfinal = maxit
else
   qfinal = q + 1
endif

  rhm1(1) = rhf(rm1)
  rhm1(2) = zhf(zm1)
  rhm1(3) = phi(phim1)

  rhm2(1) = rhf(rm1)
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
          pot_d = 0.5*(pot_it(rd,zd,phid) + pot_it(rd+1,zd,phid))
          pot_e = 0.5*(pot_it(re,ze,phie) + pot_it(re+1,ze,phie))
          
          psi_a = 0.5*(psi(ra,phia) + psi(ra-1,phia))          
          psi_b = 0.5*(psi(rb,phib) + psi(rb+1,phib))                    
          psi_c = 0.5*(psi(rc,phic) + psi(rc+1,phic))
          psi_d = 0.5*(psi(rd,phid) + psi(rd+1,phid))
          psi_e = 0.5*(psi(re,phie) + psi(re+1,phie))
          

          omsq(qfinal) = - (pot_a-pot_b)/(psi_a-psi_b)  	  

  ! and calculate the two integration constants
          c1(qfinal) = pot_b + omsq(q)*psi_b
          c2(qfinal) = pot_c + omsq(q)*psi_c
       
          rho_1d = 0.5*(rho(rd,zd,phid) + rho(rd+1,zd,phid))
          rho_c1d=rho_1d*muc1/mu1       
          h_e1d = c1(q) - pot_d - omsq(q)*psi_d
          h_c1d = h_e1d * (nc1+1)/(n1+1)*mu1/muc1                    
          cc1(qfinal) = h_c1d + pot_d+ omsq(q)*psi_d 
          
          rho_2e = 0.5*(rho(re,ze,phie) + rho(re+1,ze,phie))
          rho_c2e=rho_2e*muc2/mu2       
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
  
    
  call binary_output(c1, c2, cc1, cc2, omsq, hm1, hm2, mass1, mass2, psi, h, &
	            qfinal, initial_model_type, model_number, ra, za, phia,  &
	            rb, zb, phib, rc, zc, phic, rd, zd, phid, re, ze, phie,  &
	            rhm1, rhm2, rhom1, rhom2, xavg1, xavg2, separation, &
	            com, volume_factor, eps, hem1, hem2, rhoem1, rhoem2,     &
	            mass_c1, mass_c2, rho_1d, rho_c1d, rho_2e, rho_c2e)

  call output('density.bin', rho)

  call cpu_time(time2)

  write(13,*) 'Model: ', model_number, ' done in time: ', time2 - time1
  close(13)

end subroutine binary_scf
