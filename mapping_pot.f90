       subroutine mapping_pot(xe,qtb,x0,p0,fnon,ftb) 
         use global
         use pimc
         implicit none

         integer i,j,k,m
         real*8, intent(in)::  xe(ndof,nb)
         real*8, intent(in)::  x0(nstate,nb),p0(nstate,nb),qtb(nbath)
        
         real*8,intent(out)::  fnon(ndof,nb)
         real*8, intent(out) :: ftb(nbath)
       
         ! mapping hamiltonian variables
         real*8  hel(nstate,nstate,nb)
         real*8 hel0
          ! bead average mapping hamiltonian
         real*8 hel_ba(nstate,nstate)
         real*8  dhel(nstate,nstate,ndof,nb)
         real*8  dhel0(ndof)
         real*8 :: sum1

         fnon = 0.d0

         ! state-indepednt term for solvent
         sum1 = 0.0d0

         do i = 1,nbath
            sum1 = sum1 + (xe(1,1)/mnuc)*(ctb(i)**2/omega_b(i)**2)-ctb(i)*qtb(i)
         enddo 

         ! state dependent potential 
       
         do k = 1,nb

            ! get state dependent potential and its derivative, for each bead dimension
            call gethel(xe(:,k),hel(:,:,k))

            ! propagate mapping varibale
            call mapverlet(hel(:,:,k),x0(:,k),p0(:,k))

           ! state dependent part of the force
            call getdhel(xe(:,k),dhel(:,:,:,k))

            do m=1,ndof
               do i = 1,nstate
                  do j = 1,nstate
                     
                     fnon(m,k) = fnon(m,k) - 0.5*dhel(i,j,m,k)*(x0(i,k)*x0(j,k)+p0(i,k)*p0(j,k)) !-krondel(i,j))
                     
                  end do
               end do
            end do

            ! End the bead loop
         end do
         
         fnon(1,1) = fnon(1,1) - sum1

         !====== bath force========
         do k = 1,nbath

             ftb(k) = -mnuc*(omega_b(k)**2)*qtb(k) + ctb(k)*xe(1,1)

         enddo

         !===================== 


         return
       end subroutine mapping_pot
