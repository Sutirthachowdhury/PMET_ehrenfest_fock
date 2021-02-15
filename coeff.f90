subroutine coeff(pot,coeff1)

  use global, only : nb,nstate,ndof,betap,eigen,EVECTOR !hel1
  use pimc, only : coeff2

  implicit none

!  real*8,intent(in) :: pot_mc(:,:)
!  real*8,intent(out) :: coeff1_mc(:,:),coeff2(:,:)

  real*8, intent(in) :: pot(nstate,nstate)
  real*8, intent(out) :: coeff1(nstate,nstate)
!  real*8 ::  hel(nstate,nstate)

!  real*8, intent(out):: coeff1_mc(nstate,nstate),coeff2_mc(nstate,nstate)

   integer :: i,j

!  pot(1,1) = 0.0d0; pot(2,2) = 0.0d0; pot(1,2) = 0.0d0

  !*** v11
!  do i = 1 , ndof
!     pot(1,1) = pot(1,1) + 0.5*k1*((r(i)-r1)**2) + eps1
     !*** v22
!     pot(2,2) = pot(2,2) + 0.5*k2*((r(i)-r2)**2) + eps2
     !*** V12
!     pot(1,2) = pot(1,2) + c*exp(-alpha*((r(i)-r12)**2))
!  end do

!  pot(2,1) = pot(1,2)


  do i=1,nstate
     coeff2(i,i)=0.5*pot(i,i)/dble(nb)*dexp(-0.5*betap*pot(i,i))
  !   coeff1(i,i)=dexp(-0.5*betap*pot(i,i))
     do j=1,nstate
      !  if(j.ne.i) coeff1(i,j)=-0.5*betap*pot(i,j)*&
           !  dexp(-0.5*betap*pot(j,j))
        if(j.ne.i) coeff2(i,j)=0.5*pot(i,j)/dble(nb)*&
             dexp(-0.5*betap*pot(j,j))*(-0.5*betap*pot(j,j)+1.0d0)
     end do
  end do


!  write (*,*) "hel=",pot

!  pot(1,1) = 1.
!  pot(2,2) = 1.
!  pot(1,2)= -1.
!  pot(2,1)= -1.

  call DIAG (eigen,EVECTOR,pot)

!  write(*,*) "eigen=",eigen
!  write(*,*) "evec=",EVECTOR (1,1), EVECTOR (1,2)
!  write(*,*) "evec=",EVECTOR (2,1), EVECTOR (2,2)


  coeff1 = 0.

  do i=1,nstate
     do j=1,nstate
        coeff1(i,j) = sum(EVECTOR(i,:)*dexp(-betap*eigen(:))*EVECTOR(j,:))
     enddo
  enddo

!  do i=1,nstate
!     write(*,*) (coeff1 (i,j),j=1,nstate)
!  end do


!========================================================================================!
! Ananth-Miller approximated version (JCP2010) based on Chandler's short time approximation
!  do i=1,nstate
!     coeff2(i,i)=0.5*pot(i,i)/dble(nb)*dexp(-0.5*betap*pot(i,i))
!     coeff1(i,i)=dexp(-0.5*betap*pot(i,i))
!     do j=1,nstate
!        if(j.ne.i) coeff1(i,j)=-0.5*betap*pot(i,j)*&
!             dexp(-0.5*betap*pot(j,j))
!        if(j.ne.i) coeff2(i,j)=0.5*pot(i,j)/dble(nb)*&
!             dexp(-0.5*betap*pot(j,j))*(-0.5*betap*pot(j,j)+1.0d0)
!     end do
!  end do
!=========================================================================================!

  return

end subroutine coeff
