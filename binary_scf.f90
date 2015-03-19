subroutine binary_scf(model_number, initial_model_type, ra, rb, rc, rd, re, rhom1, rhom2, frac, pin, eps, &
                     qfinal)
implicit none
include 'runscf.h'
!include 'mpif.h'
!**************************************************************************************************
!
!  subroutine arguments
!

integer, intent(in) :: model_number
integer, intent(in) :: initial_model_type
integer, intent(in) :: ra, rb, rc, rd, re
real, intent(in) :: rhom1
real, intent(in) :: rhom2
real, intent(in) :: frac
real, intent(in) :: pin
real, intent(in) :: eps
integer, intent(out) :: qfinal

!
!**************************************************************************************************
!
!   global varaibles
!

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd) :: r, rhf, rinv, rhfinv
real, dimension(numz_dd) :: zhf
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
!**************************************************************************************************
!
!    local variables
!

real, dimension(numr_dd, numz_dd, numphi) :: h, pot_it, pot_old, temp

real, dimension(numr_dd, numphi) :: psi

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
logical :: temp_processor

integer :: ra_pe, rb_pe, rc_pe

logical :: i_have_ra, i_have_rb, i_have_rc

integer :: ra_local_index, rb_local_index, rc_local_index

real :: temp_hm1, temp_hm2

real, dimension(3) :: rhm1, rhm2, temp_rhm1, temp_rhm2

integer :: phi1, phi2, phi3, phi4

integer :: I, J, K, L, Q

integer :: ierror

  real :: pot_a, pot_b, pot_c, pot_d, pot_e
  real :: psi_a, psi_b, psi_c, psi_d, psi_e


  real :: rho_c1d, rho_1d, rho_c2e, rho_2e
  real :: h_c1d, h_e1d, h_c2e, h_e2e 
  real :: rhoem1, rhoem2, norm1, norm2
  real :: rem1, zem1, phiem1, rem2, zem2, phiem2
  real :: rm1, zm1, phim1, rm2, zm2, phim2
  integer :: rmax
   real :: gammac1, gammac2, gammae1,gammae2
   real :: kappac1,kappac2, kappae1,kappae2

!integer, dimension(MPI_STATUS_SIZE) :: istatus 

!
!**************************************************************************************************

call cpu_time(time1)

qfinal = 1

phi1 = int(numphi / 4.0) - 1
phi2 = int(numphi / 4.0) + 1
phi3 = int(3.0 * numphi / 4.0) - 1
phi4 = int(3.0 * numphi / 4.0) + 1

   gammae1 = 1.0 + 1.0/n1
   gammae2 = 1.0 + 1.0/n2
   gammac1 = 1.0 + 1.0/nc1
   gammac2 = 1.0 + 1.0/nc2

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

i_have_ra = .false.
ra_local_index = 0
ra_pe = 0

if ( iam_on_bottom ) then
   do I = rlwb, rupb
      if ( abs(rhf_g(ra) - rhf(I)) < 0.1 * dr ) then
         i_have_ra = .true.
         ra_local_index = I
      endif
   enddo
endif


i_have_rb = .false.
if ( iam_on_bottom ) then
   do I = rlwb, rupb
      if ( abs(rhf_g(rb) - rhf(I)) < 0.1 * dr ) then
         i_have_rb = .true.
         rb_local_index = I
      endif
   enddo
endif



i_have_rc = .false.
if ( iam_on_bottom ) then
   do I = rlwb, rupb
         if ( abs(rhf_g(rc) - rhf(I)) < 0.1 * dr ) then
          i_have_rc = .true.
          rc_local_index = I
        endif
   enddo
endif


volume_factor = 2.0 * dr * dz * dphi
  rmax = numr - 8
!gamma = 1.0 + 1.0 / pin
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

if ( iam_root ) then
   open(unit=13,file='iteration_log',form='formatted',status='unknown',position='append')
   write(13,*) iam, 1, mass1(1), mass2(1), xavg1, xavg2, com, separation
endif

print*, 'maxit = ', maxit

do Q = 2, maxit-1                                   ! START OF THE ITERATION CYCLE
print*, "================"
print*, "iteration number = ", Q
   ! solve the Poisson equation for the current density field
   if ( Q == 2 ) then
      call potential_solver(0)
      pot_old = pot
   else
      call potential_solver(1)
   endif
   pot_it = (1.0 - frac) * pot + frac * pot_old
   pot_old = pot

   ! compute the form factor for the centrifugal potential
   do K = philwb, phiupb
      do I = rlwb-1, rupb+1
         psi(I,K) = - 0.5 * ( (rhf(I)*cosine(K) - com)**2 + rhf(I)*rhf(I)*sine(K)*sine(K) )
      enddo
   enddo

   ! calculate the angular frequency
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
         
!print*, pot_a, pot_b, pot_c, pot_d, pot_e
!print*, psi_a, psi_b, psi_c, psi_d, psi_e 

          omsq(q) = - (pot_a-pot_b)/(psi_a-psi_b)
          
          c1(q) = pot_b + omsq(q)*psi_b
          c2(q) = pot_c + omsq(q)*psi_c
       
          rho_1d = 0.5*(rho(rd,zd,phid) + rho(rd+1,zd,phid))
          rho_c1d=rho_1d*muc1/mu1       
          h_e1d = c1(q) - pot_d - omsq(q)*psi_d
          h_c1d = h_e1d * (nc1+1)/(n1+1)*mu1/muc1                    
          cc1(q) = h_c1d + pot_d+ omsq(q)*psi_d 

print*, "rho_1d", rho_1d, "rho_c1d", rho_c1d       

          rho_2e = 0.5*(rho(re,ze,phie) + rho(re+1,ze,phie))
          rho_c2e=rho_2e*muc2/mu2       
          h_e2e = c2(q) - pot_e - omsq(q)*psi_e
          h_c2e = h_e2e * (nc2+1)/(n2+1)*mu2/muc2                    
          cc2(q) = h_c2e + pot_e+ omsq(q)*psi_e       
print*, "rho_2e", rho_2e, "rho_c2e", rho_c2e

!print*, "rhos", rho_c1d, rho_1d, rho_c2e, rho_2e
!print*,  omsq(q) 
!print*,  "cs",c1(q), c2(q),cc1(q),cc2(q)

 
   ! now compute the new  enthalpy field from the potential and SCF constants
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

!print*, "post loop"
!print*, "h's",hem1(q), hm1(q)
!print*,"norms",norm1, norm2 


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
            if ( rhf(I) <= rhf_g(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi2, phi3
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf_g(rc) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo
   do K = phi4, phiupb
      do J = zlwb, zupb
         do I = rlwb, rupb
            if ( rhf(I) <= rhf_g(rb) ) then
               rho(I,J,K) = 0.0
            endif
         enddo
      enddo
   enddo

   ! impose the equatorial boundary condition
   if ( iam_on_bottom ) then
      do K = philwb, phiupb
         do I = rlwb, rupb
            rho(I,zlwb-1,K) = rho(I,zlwb,K)
         enddo
      enddo
   endif

   ! impose the axial boundary condition
   if ( iam_on_axis ) then
      rho(rlwb-1,:,:) = cshift(rho(rlwb,:,:),dim=2,shift=numphi/2)
   endif

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
   call compute_virial_error(psi, h, sqrt(omsq(Q)), pin, volume_factor, virial_error1, &
                             virial_error2, virial_error)

print*, "omsq = ", cnvgom
print*, "c1 = ",  cnvgc1,  "c2 = ",  cnvgc2
print*,  "h1 = ",  cnvgh1,  "h2 = ",  cnvgh2
  kappac1 = rhom1*hm1(q)/(nc1+1.0)/rhom1**(gammac1)
  kappac2 = rhom2*hm2(q)/(nc2+1.0)/rhom2**(gammac2)
  kappae1 = kappac1*rho_c1d**gammac1/rho_1d**gammae1
  kappae2 = kappac2*rho_c2e**gammac2/rho_2e**gammae2

print*, "hm1(q) = ",hm1(q),"hm2(q) = ",hm2(q)
print*, "rmom1 = ", rhom1, "rhom2 = ", rhom2
!print*, "",,"",

 print*, "kappac1= ", kappac1,  "kappae1", kappae1
 print*, "kappac2= ", kappac2,  "kappae2", kappae2

   if ( iam_root ) then
      write(13,*) Q, mass1(Q), mass2(Q), xavg1, xavg2, com, omsq(Q), c1(Q), c2(Q), &
                  hm1(Q), hm2(Q), cnvgom, cnvgc1, cnvgc2, cnvgh1, cnvgh2, virial_error1, &
                  virial_error2, virial_error
   endif

   if ( cnvgom < eps .and. cnvgc1 < eps .and. cnvgc2 < eps .and. &
        cnvgh1 < eps .and. cnvgh2 < eps ) then
      exit
   endif

!   if ( virial_error > virial_error_prev .and. Q > 10  ) then
!      exit
!   endif

enddo                                               ! END OF THE ITERATION CYCLE

if ( q >= maxit ) then
   qfinal = maxit
else
   qfinal = q + 1
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
  
if ( iam_root ) then
   write(13,*) iam, 'Model: ', model_number, ' done in time: ', time2 - time1
endif

  call binary_output(c1, c2, cc1, cc2, omsq, hm1, hm2, mass1, mass2, psi, h, &
	            qfinal, initial_model_type, model_number, ra, za, phia,  &
	            rb, zb, phib, rc, zc, phic, rd, zd, phid, re, ze, phie,  &
	            rhm1, rhm2, rhom1, rhom2, xavg1, xavg2, separation, &
	            com, volume_factor, eps, hem1, hem2, rhoem1, rhoem2,     &
	            mass_c1, mass_c2, rho_1d, rho_c1d, rho_2e, rho_c2e)


  call ancient_output(c1, c2, omsq, hm1, hm2, mass1, mass2, psi, h, qfinal,     &
                   initial_model_type, model_number, ra, za, phia, rb, zb,   &
                   phib, rc, zc, phic, rhm1, rhm2, 1.5, rhom1, rhom2, xavg1, &
                   xavg2, separation, com, volume_factor, eps)



call cpu_time(time1)

call output(1000, 'density.bin', rho)

call cpu_time(time2)

if ( iam_root ) then
   write(13,*) iam, 'Model: ', model_number, ' disk I/O done in time: ', time2 - time1
   close(13)
endif

end subroutine binary_scf
