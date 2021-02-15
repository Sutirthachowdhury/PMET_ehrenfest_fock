subroutine init_nuclear_sampling(xe_mc)
 

  use global
  use pimc
  !use globran

  implicit none
 
  real*8, intent(out) :: xe_mc(ndof,nb)
  real*8 :: gran
  real*8 :: drp0,drp,eta
  integer :: i,k,l
  
  
  do i=1,ndof

     drp = sqrt(1.0/(beta*mnuc*omega(i)**2))

     xe_mc(i,1) = gran()*drp

  end do

  
 end subroutine init_nuclear_sampling
  

  
