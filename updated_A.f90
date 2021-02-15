subroutine updated_A(xe_new,ve_new,x0_new,p0_new,estimator)
  
use global
use pimc

implicit none

real*8 :: hel(nstate,nstate),pop_init_nu,pop_init_de

real*8, intent(in) :: xe_new(ndof,nb),ve_new(ndof,nb),x0_new(nstate,nb),p0_new(nstate,nb)

real*8, intent(out) :: estimator

integer ::  i,j,k,l,iup,idn

estimator = 0.0d0

do i = 1,nb
   
   call gethel (xe_new(:,i),hel)
   call coeff(hel)



   iup=i+1
   if(iup.gt.nb) iup=1
   idn=i-1
   if(idn.lt.1) idn=nb

   pop_init_nu = 0.
   pop_init_de = 0.

   do l = 1,nstate
      pop_init_nu = pop_init_nu + (p0_new(1,idn)*coeff1(1,l)*x0_new(l,i))
      do k = 1,nstate
         pop_init_de = pop_init_de + (p0_new(l,idn)*coeff1(l,k)*x0_new(k,i))
      enddo
   enddo

   estimator = estimator + (pop_init_nu/pop_init_de)


enddo

estimator = estimator/real(nb)
!write(*,*) "estimator=", estimator

end subroutine updated_A
