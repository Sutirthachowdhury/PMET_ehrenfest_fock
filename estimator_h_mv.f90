subroutine estimator_h_mv(x0,p0,theta)

use global,only : nb,nstate,ndof,beta,betap,eye

implicit none

real*8, intent(in) :: x0(nstate,nb),p0(nstate,nb)
complex*16 :: gamma1(nstate,nstate)
complex*16,intent(out) :: theta
integer :: i,iup,idn,j,k,l
  
real*8 :: sums
real*8 :: identity(nstate,nstate)
complex*16 :: weight(nstate,nstate)

!definition of 3*3 identity matrix
identity(1,1)=1.
identity(1,2)=0.
identity(2,1)=identity(1,2)
identity(2,2)=identity(1,1)
!----------------------------------------
!definition of initial gamma matrix--------
gamma1(1,1) = cmplx(1.,0.)
gamma1(1,2) = cmplx(0.,0.)
gamma1(2,1) = gamma1(1,2)
gamma1(2,2) = gamma1(1,1)
!----------------------------------------
do i=1,nb

   iup=i+1
   if(iup.gt.nb) iup=1
   idn=i-1
   if(idn.lt.1) idn=nb

   sums = sum((x0(:,i)**2 + p0(:,i)**2))

   !initialization of weight
   weight(:,:)=cmplx(0.,0.)

   do j=1,nstate
      do k=1,nstate
         !for prefactor
         weight(j,k) = dexp(-sums)*(((x0(j,i)+eye*p0(j,i))*(x0(k,i)-eye*p0(k,i)))-0.5*identity(j,k))
      enddo
   enddo
   !full prefactor
   gamma1 = MATMUL(gamma1,weight)
  ! write(*,*) gamma1(1,1),gamma1(2,2)
enddo

!write(*,*) gamma1(1,1),gamma1(2,2)
!=================== taking trace===============================
theta=cmplx(0.,0.)

do l = 1,nstate
   theta = theta + gamma1(l,l)
enddo
!=============================================================

end subroutine estimator_h_mv
