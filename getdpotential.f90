subroutine getdpotential(qs,dhel)

  use global
  use pimc
  implicit none

  integer i,j,l

  real*8, intent(in) :: qs(ndof)
  
  real*8, intent(out) :: dhel(nstate,nstate,ndof)

  dhel(1,1,1)=0.0d0 ; dhel(2,2,1)=0.0d0 ; dhel(1,2,1)=0.0d0 ; dhel(2,1,1)=0.0d0
!============== derivative of model 1 ananth JPCL ===============================
  !**derivative of V11
  do i = 1,ndof
     if(qs(i).lt.0) then
        dhel(1,1,i) = 0.016*dexp(1.6*qs(i))
     else
        dhel(1,1,i) = 0.016*dexp(-1.6*qs(i))
     endif

  !**derivative of V12
     dhel(1,2,i) = -0.01*qs(i)*dexp(-qs(i)**2)
  enddo

  dhel(2,2,1) = -dhel(1,1,1)

  dhel(2,1,1) = dhel(1,2,1)

!==============================================================================

 !===============derivative of model 2 ananth JPCL =============================

!  dhel(1,1,1) = 0.0
!  do i = 1,ndof
!     dhel(2,2,i) = 0.056*qs(i)*dexp(-0.28*qs(i)**2)
!     dhel(1,2,i) = -0.0018*qs(i)*dexp(-0.06*qs(i)**2)
!  enddo
!  dhel(2,1,1) = dhel(1,2,1)
  !=============================================================================

!==================derivative of model 3 ananth JPCl===========================
  !dhel(1,1,1) = 0.0
  !dhel(2,2,1) = 0.0
  
 ! do i = 1,ndof
  !   if(qs(i).lt.0) then
   !     dhel(1,2,i) = 0.0018*dexp(0.9*qs(i))
    ! else
     !   dhel(1,2,i) = 0.0018*dexp(-0.9*qs(i))
     !endif
  !enddo

  !dhel(2,1,1) = dhel(1,2,1)
!================================================================================ 
  return
end subroutine
