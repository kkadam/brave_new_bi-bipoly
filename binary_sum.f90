subroutine binary_sum(temp, ret1, ret2)
implicit none
include 'runscf.h'

real, intent(in), dimension(numr_dd,numz_dd,numphi) :: temp
real, intent(out) :: ret1, ret2
integer :: phi1, phi2, phi3, phi4
integer :: I, J, K

phi1 = int(numphi / 4.0) - 1
phi2 = int(numphi / 4.0) + 1
phi3 = int(3.0 * numphi / 4.0) - 1
phi4 = int(3.0 * numphi / 4.0) + 1

ret1 = 0.0
ret2 = 0.0

do K = philwb, phi1
   do J = zlwb, zupb
      do I = rlwb, rupb
         ret1 = ret1 + temp(I,J,K)
      enddo
   enddo
enddo
do K = phi2, phi3
   do J = zlwb, zupb
      do  I = rlwb, rupb
         ret2 = ret2 + temp(I,J,K)
      enddo
   enddo
enddo
do K = phi4, phiupb
   do J = zlwb, zupb
      do I = rlwb, rupb
         ret1 = ret1 + temp(I,J,K)
      enddo
   enddo
enddo

return
end subroutine binary_sum
