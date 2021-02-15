subroutine parameter_input
      use global 
      use pimc
      implicit none
 
      integer i
      real*8 w0,wm
      !pi = dacos(-1.d0)
      !-------------------------------------------------------
      open(10,file='params.in',status='old')
      read(10,*) nb
      read(10,*) nstep
      read(10,*) ntraj
      read(10,*) beta
      read(10,*) masse
      read(10,*) xi
      read(10,*) wc
      read(10,*) mnuc
      read(10,*) omega_c
      read(10,*) delta
      read(10,*) f0
      read(10,*) shift
      read(10,*) bias
      read(10,*) dt
      read(10,*) ourseed
      close(10)
      !---------------------------------------------------------


      allocate(delta_r0(nbox),step(nstate+1),pop(nstate,nbox),population(nbox),inv(nstate,nstate),coeff2(nstate,nstate))
      
     ! allocate(eigen(nstate),EVECTOR(nstate,nstate),hel1(nstate,nstate))
      allocate(eigen(nstate),EVECTOR(nstate,nstate))
      
      betap = beta/dble(nb)
      print *, 'betap', betap

 
      open(10,file='steps.in',status='old')
      read(10,*) step(1), step(2), step(3)
      close(10)
 
      ! hard coded, need to be modified later, with setup subroutine      

      omega(:) = 3.50690243E-4



!========== spectral density for continuous bath ===================

      wm=3.0*wc ! bath maximum freq                                   

      w0=wc*(1.-exp(-2.0))/real(nbath)

      do i=1,nbath
         omega_b(i) = -wc*log(1.0-i*w0/wc)

         ctb(i) =  sqrt(xi*w0*mnuc)*omega_b(i)
      
      end do
!===================================================================


      ! for mapping hamiltonian                                                                                              

       krondel=0.

       do i=1,nstate
          krondel(i,i)=1.
       end do

end subroutine parameter_input
