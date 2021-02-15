subroutine plot_potential
  use global,only : nstate,ndof 
  implicit none
  real*8 :: xe_mc(ndof)
  real*8 :: hel(nstate,nstate),hel0
  real*8 :: dhel(nstate,nstate,ndof)
  integer :: i
 
  do i=1,200
     xe_mc= -0.15 + i*0.0015d0
     call gethel(xe_mc,hel)
     write(12,*) xe_mc, hel(1,1), hel(2,2)
     write(13,*) xe_mc, hel(1,2)
  end do
  return
end subroutine plot_potential
