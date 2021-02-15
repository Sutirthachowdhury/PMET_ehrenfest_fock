subroutine init_mapping_sampling(x0_mc,p0_mc)

use global
!use globran
implicit none

real*8,intent(out) :: x0_mc(nstate,nb),p0_mc(nstate,nb)
real*8 :: d_mom,d_pos,gran
real*8 :: Rfocus(nstate),rand,theta
integer :: i,j,k,init_state 


!-----------------------
do j = 1,nstate
   x0_mc(j,1) = 0.0
   p0_mc(j,1) = 0.0
enddo
!-----------------------

!--- initial states-----
x0_mc(1,1) = 1.0
p0_mc(1,1) = 1.0
!----------------------
end subroutine init_mapping_sampling

    
