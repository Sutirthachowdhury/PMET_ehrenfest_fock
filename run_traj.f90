       Subroutine run_traj(xe,ve,qtb,vtb,x0,p0,fnon,ftb) 
         use global
         use pimc
         !use globran
         !use mpi_lib


         implicit none

         integer i,j,k,m
 
         real*8, intent(inout):: xe(ndof,nb),ve(ndof,nb),x0(nstate,nb),p0(nstate,nb),fnon(ndof,nb)
         real*8, intent(inout) :: qtb(nbath), vtb(nbath) , ftb(nbath)
         
         
         
         real*8 :: f1(ndof,nb),f2(ndof,nb),aa(ndof,nb)

         real*8 invmasse,sigma,gama,gran,r1,r2,taul

         integer istep

         real*8  dt2

         !!!!!!!!!!!!!! Initialize correlation function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         taul = 19139.8357

         dt2  = dt*0.5d0

         invmasse = 1.d0 / mnuc

        !-------------------------------------------
         ve(1,1) = ve(1,1) + invmasse * fnon(1,1)*dt2
            
         do k = 1,nbath
           vtb(k) = vtb(k) + invmasse * ftb(k)*dt2
         enddo

        !--------------------------------------------------

         ! Evolve the ring polymer according to its internal potential...
         call mapping_pot(xe,qtb,x0,p0,fnon,ftb) 
         
        
         xe(1,1) = xe(1,1) + ve(1,1)*dt !+ 0.5 * invmasse * fnon*dt**2
         
         do k = 1,nbath
          qtb(k) = qtb(k) + vtb(k)*dt  !+ 0.5*(ftb(k)/mnuc)*dt**2
         enddo
     

         call mapping_pot(xe,qtb,x0,p0,fnon,ftb)

          !-----------------------------------------------
          
         do k = 1,nbath
           vtb(k) = vtb(k) + invmasse * ftb(k)*dt2
         enddo
        
         
         ve(1,1) = ve(1,1) + invmasse * fnon(1,1)*dt2
         


 
222      format(9(f13.6,2x))
       end subroutine run_traj
