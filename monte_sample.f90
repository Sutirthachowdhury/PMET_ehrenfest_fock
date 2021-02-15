subroutine monte_sample(xe_mc)

  use pimc, only : delta_r0,step,nmc,nbox,step,pop,ef_eta,&
       avg_energy,mnuc
  use global,only : nb,nstate,ndof,beta,betap,pi

  implicit none

      !index for which eta value
      integer :: itype,ibead,ifile,imove,nfile
      integer :: int_nmc,nstore,stat_count
      integer :: i,j,iup,idn, k, int_mc

      real*8 :: iacc_x,iatt_x,iacc_s,iacc_s1,iatt_s,iatt_s1, rnum,step_mc
      real*8 :: wt,wt1,wtt,eterm,nterm
      real*8 :: avg_part,std_dev_en,std_dev_part,ecorr,pcorr,error,sums
      real*8 :: avg_sign,dr
      real*8 :: plus_energy, minus_energy, plus_part, minus_part
 
      real*8, allocatable:: xe_dn(:),xe_up(:),xe_new(:,:),x0_new(:,:),p0_new(:,:)
      real*8, intent(inout):: xe_mc(ndof,nb)
      
      !real*8, intent(out):: pref !inv1(nstate,nstate) 
      
      complex*16 :: gamma
      complex*16 :: gamma1
      
      !real*8,allocatable  :: state1(:),state2(:),statenew(:),P1(:),P2(:),&
                           !Pnew(:),store1(:,:),store2(:,:,:),store_en(:),&
                           !store_part(:),avg1(:),avg2(:,:),std_dev1(:),std_dev2(:,:),pn(:)

       real*8, allocatable :: store1(:,:),store2(:,:,:),store_en(:), store_part(:)
       real*8, allocatable :: avg1(:),avg2(:,:),std_dev1(:),std_dev2(:,:),pn(:)


       !random number used in MC                                                                                                            
      real*8 :: rand

      !timing the program run                                                                                                                                                           
      real :: etime, elapsed(2), total

      !no: of blocks for data averaging                                                                                                                                                
      nfile=1
      nstore=200


      allocate(xe_dn(ndof),xe_up(ndof),xe_new(ndof,nb),x0_new(nstate,nb),p0_new(nstate,nb))

      !allocate(state1(nstate),state2(nstate),statenew(nstate))

      !allocate(P1(nstate),P2(nstate),Pnew(nstate))

      allocate(store1(nfile,nbox),store2(nfile,nstate,nbox),store_en(nfile),store_part(nfile))

      allocate(avg1(nbox),avg2(nstate,nbox),std_dev1(nbox),std_dev2(nstate,nbox),pn(nstate))

      !=========MONTE CARLO INTEGRATION RUN==========                                                                                                                                   

      !initialize all zeros                                                                                                                                                             

      avg_sign=0.d0
      iatt_x=0.d0
      iacc_x=0.d0
      iatt_s=0.d0
      iatt_s1=0.d0
      iacc_s=0.d0
      iacc_s1=0.d0

      avg_energy=0.d0
      avg_part=0.d0
      delta_r0=0.d0
      pop=0.d0
      int_nmc=1
      stat_count=0
      plus_part=0.d0;plus_energy=0.d0
      minus_part=0.d0;minus_energy=0.d0

      !ifile loop key                                                                                                                                                                  
         step_mc = 1.d0
!         int_nmc = 1
         !stepmc loop key                                                                                                                                                               
         do int_mc=1,nmc
!            write(*,*) "xe_mc",xe_mc 
            !loop over all beads in each mc step                                                                              
!            write(500,222) xe_mc(1,1), x0_mc(1,1),x0_mc(2,1),p0_mc(1,1),p0_mc(2,1)                                                         
 
            do i=1,nb

               !randomly pick a bead to move                                                                                                     
               call random_number(rand)
               ibead = int(rand*nb)+1

               iup=ibead+1
               if(iup.gt.nb) iup=1
               idn=ibead-1
               if(idn.lt.1) idn=nb
!               giev the currnt xe_mc to xe_new
               xe_new = xe_mc

               do j=1,ndof
                  call random_number(rand)
                  rnum=rand 
                  xe_new(j,ibead) = xe_mc(j,ibead)+(rnum-0.5d0)*step(1)
                  xe_dn(j)=xe_mc(j,idn);xe_up(j)=xe_mc(j,iup)
!                  write(*,*) "xnew=",xe_new,"x1=",xe_dn,"x2=",xe_up
               end do

!               write(*,*) "xe_mc",xe_mc 
               call sampling_nuc(xe_dn,xe_new(:,ibead),xe_up,wtt)

!                  iatt_x=iatt_x+1.0
                 
          !write(*,*) "imove=",imove,"xnew=",xnew,"statenew=",statenew

                  !sampling call key
                  !checking against original weight

 
                  xe_dn(:)=xe_mc(:,idn);xe_up(:)=xe_mc(:,iup)

                  call sampling_nuc(xe_dn,xe_mc(:,ibead),xe_up,wt1)
                  
                  !accepting/rejecting step

                  wt=wtt/wt1
                  call random_number(rand)
                  rnum=rand

                  !write(*,*) 'wt=' , wt

                  if(rnum.lt.wt) then
                     xe_mc(:,ibead)=xe_new(:,ibead)
!                     iacc_x=iacc_x+1.d0
                  end if


               !end loop over nb beads
               end do

              !end loop over nmc                                                                              
            end do

     !print *,'acceptance rate',iacc_x/iatt_x,iacc_s/iatt_s,iacc_s1/iatt_s1
    
!        open(20,file='final_config.dat')
!         write(20,*) xe_mc
!         write(20,*) x0_mc
!         write(20,*) p0_mc
!        close(20)

222 format(6(f13.6,2x))                                                                                                                   
!Return                                                                                                                                    
!end subroutine monte_sample                                                                                                               
end subroutine monte_sample
!----------------------------------------------------------------------------------------------------------------------                    

!===============================================================                                                                           
!                   NUCLEAR SAMPLING                                                                                                       
!===============================================================                                                                           
subroutine sampling_nuc(x1,x2,x3,wtt)

          use pimc, only : inv,mnuc !,eye
          use global,only : nstate,nb,ndof,beta,betap,masse,omega,eye

          implicit none

          real*8,intent(in) :: x1(ndof),x2(ndof),x3(ndof)
          real*8,intent(out) :: wtt
   
          real*8 :: hel(nstate,nstate)

          integer :: i,j,k
          real*8 :: const,pref,term,sums1,sums2

          wtt=1.d0
          const=0.5d0*mnuc/betap
          !ke part of the nuclei(ring polymer part only for more than 1 bead)                                                                                                           
          wtt=wtt*dexp(-const*sum((x1(:)-x2(:))**2+(x2(:)-x3(:))**2))
         ! state independent part of the potential
          ! XXXXX modify here 
!          call getu0(x2,u0)

          ! for spin-boson 
          wtt=wtt*dexp(-betap*sum(0.5*omega(:)**2*x2(:)**2))

          !for Ananth JCP(2013)
          !wtt=wtt*dexp(-betap*sum(0.5*k1*x2(:)**2))
!          wtt=wtt*dexp(-betap*u0)
          !call potential(x2,hel)
         
         ! call gethel (x2,hel)
          !call coeff(hel)
         
          !inv1=coeff1                                                                                                                     

          !gaussian term form e-ic variables and prefactor term                                                                            
          !sums1=0.0d0
          !sums2=0.0d0
          !do j=1,nstate
           !  do k=1,nstate
                !gaussian part                                                                                                             
            !    sums1=sums1+state2(j)*coeff1(j,k)*P2(k)
             !   sums2=sums2+P1(j)*coeff1(j,k)*state2(k)
             !end do
          !end do
          !wtt=wtt*dabs(sums1*sums2)
         
end subroutine sampling_nuc
!---------------------------------------------------------------------------------------------------------------------------------         



