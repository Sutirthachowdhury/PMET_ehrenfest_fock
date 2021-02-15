subroutine potential(r,hel)

  !use pimc, only : k1, k2, r1, r2, eps1, eps2, c, gamma, r12
  use pimc
  use global,only : nstate,ndof,nb
 
  implicit none

  real*8, intent(in) :: r(ndof)
  real*8, intent(out) :: hel(nstate,nstate)
 
  integer :: i,j

  hel(1,1) = 0.0d0 ; hel(2,2) = 0.0d0 ; hel(1,2) = 0.0d0 ; hel(2,1) = 0.0d0

!===============Ananth JPCL model1 =========================== 
     !*** v11
  do i = 1 , ndof
     if(r(i).lt.0) then
       hel(1,1) =  (-0.01*(1-dexp(1.6*r(i))))+0.01     
    else
       hel(1,1) =  (0.01*(1-dexp(-1.6*r(i))))+0.01
    endif
     !***v12
    hel(1,2) =  (0.005*dexp(-r(i)**2))
 end do

 hel(2,2) = -hel(1,1)+0.02
 hel(2,1) = hel(1,2)


!===============================================================

!==============Ananth JPCL model 2 =============================
  !*** v11
!  hel(1,1) = 0.0
  !** v22
!  do i = 1 , ndof
!     hel(2,2) = -0.1*dexp(-0.28*r(i)**2) + 0.05
  !** v12
!     hel(1,2) = 0.015*dexp(-0.06*r(i)**2)
!  end do 
  
! hel(2,1) = hel(1,2)
!================================================================

!=========== Ananth JPCL model 3 ==============================
 !*** v11
! hel(1,1) = 0.0
 !*** v22
 !hel(2,2) = 0.0012
 !*** v12
 !do i = 1 , ndof
  !   if(r(i).lt.0) then
   !     hel(1,2) =  0.002*dexp(0.9*r(i))
    ! else
     !   hel(1,2) =  0.002*(2-dexp(-0.9*r(i)))
     !endif
 ! enddo

 ! hel(2,1) = hel(1,2)
!====================================================================== 
  return

end subroutine potential
