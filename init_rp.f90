!----------------------------------------------------------------------------------------------------------------------
subroutine init_rp(xe_mc)

  use global,only : nb,nstate,ndof,omega,betap
  use pimc, only : mnuc

  implicit none
  real*8 :: rand,rnum,gran,drp
 
  real*8, intent(out) :: xe_mc(ndof,nb)

  integer :: i,j
 
 

  do i=1,ndof
     drp = 1.0/sqrt(mnuc*omega(i)**2*betap)
     
     do j=1,nb
        xe_mc (i,j) = gran()*drp + 7.0
     end do
  end do


end subroutine init_rp
!----------------------------------------------------------------------------------------------------------------------
