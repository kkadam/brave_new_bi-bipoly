subroutine binary_sum(temp, ret1, ret2)
implicit none
include 'runscf.h'

real, intent(in), dimension(numr,numz,numphi) :: temp
real, intent(out) :: ret1, ret2
integer :: I, J, K


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
