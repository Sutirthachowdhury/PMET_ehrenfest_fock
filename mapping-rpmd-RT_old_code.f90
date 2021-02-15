program rpmd
     use global
     use pimc,only:nbox,population,mnuc 

     implicit none

     
     integer :: istep,iup,idn     
     integer,allocatable :: seed(:)
     integer seed_dimension,time
   
     integer i,j,k,l,m,n,ncount
     real*8, allocatable:: xe(:,:),ve(:,:),fnon(:,:)
     real*8, allocatable:: xe_mc(:,:),x0_mc(:,:),p0_mc(:,:)
     real*8, allocatable:: x0(:,:),p0(:,:)
     real*8, allocatable:: pop_wigner(:,:) 
     real*8, allocatable :: xe_avg(:)
     real*8, allocatable :: ptb(:) , qtb(:) , ftb(:) , dhtb(:) ,vtb(:) 
     
     real*8 :: dp,eta
     real*8 :: gran,rand
     real*8 :: sigmaQtb, sigmaPtb
     real*8 :: hel(nstate,nstate)
     

     ! IO for MC input parameters
     call parameter_input
     
            
     ! Initialize the RNG seed
     !call random_seed(size=seed_dimension)
     
     !allocate(seed(seed_dimension))
     
     !do i=1,seed_dimension
     !   seed(i)=time()+3*i-1
     !enddo
     !call random_seed(put=seed)


     CALL RANDOM_SEED(size=seed_dimension)
     ALLOCATE (seed(seed_dimension))
     do i=1,seed_dimension
        seed(i) =time()+3*i-1 + ourseed
     end do
     CALL RANDOM_SEED(PUT=seed)


     allocate(xe_mc(ndof,nb),x0_mc(nstate,nb),p0_mc(nstate,nb))

     allocate(xe(ndof,nb),ve(ndof,nb),fnon(ndof,nb)) 
     allocate(x0(nstate,nb),p0(nstate,nb))
     allocate(ptb(nbath),qtb(nbath),ftb(nbath),dhtb(nbath),vtb(nbath))
     
     allocate(pop_wigner(nstate,nstep)) 
     allocate(xe_avg(nstep))

     ! initialize all varialbes to 0 /or some value
     !1call init_rp(xe_mc)

     call plot_potential
         
          
     pop_wigner = 0.
         
     xe_avg = 0.0d0


     do j=1, ntraj
            
        call init_nuclear_sampling(xe_mc)

        call init_mapping_sampling(x0_mc,p0_mc)
                   
        xe=xe_mc
        x0=x0_mc
        p0=p0_mc

        !write(3000,222) real(j),(xe(1,k),k=1,nb)
        
        !-------------velocity distribution--------------------------------     
     
       do k=1,ndof

           dp=sqrt(mnuc/beta) 

           ve(k,1) = gran()*dp

           ve(k,1) = ve(k,1)/mnuc

        end do
      !================================================

      do k=1,nbath
         sigmaQtb = sqrt(1.0/(beta*mnuc*omega_b(k)**2))
         sigmaPtb = sqrt(mnuc/beta)
         qtb(k) = gran()*sigmaQtb + ((ctb(k)*xe(1,1))/(mnuc*omega_b(k)**2))  
         ptb(k) = gran() * sigmaPtb
      end do  

      vtb(:) = ptb(:)/mnuc

         !do k = 1,nstate
         !   write(500,*) k,x0(k,1)
         !enddo

        !---------------- main dynamics propagation loop----------------------------

        call mapping_pot(xe,qtb,x0,p0,fnon,ftb)
             
        do istep=1,nstep
          
           call run_traj(xe,ve,qtb,vtb,x0,p0,fnon,ftb)

           !write(1000+j,222) real(istep*dt),(x0(k,1),k=1,nstate),(p0(k,1),k=1,nstate)
           !write(2000+j,222) real(istep*dt),(0.5*x0(k,1)**2+0.5*p0(k,1)**2,k=1,nstate)

  !         xe_avg(istep) = xe_avg(istep) + xe(1,1)

 !          write(100+j,222) real(istep),xe(1,1)

           !------------------- populations (diabatic)----------------------------------------- 
           
           do k=1,nstate
              do i=1,nb
                 pop_wigner(k,istep) =  pop_wigner(k,istep) &
                      + ((0.5*(x0(k,i)**2+p0(k,i)**2))/real(nb)) 
              enddo
           enddo
           
           !end of time step loop
        enddo
         ! end of trajectory loop
     end do
         
     pop_wigner = pop_wigner/real(ntraj)
!     xe_avg = xe_avg/real(ntraj)
                  
     do ncount=1,nstep
        write(80,222) ncount*dt,(pop_wigner(k,ncount),k=1,nstate) !,xe_avg(ncount)
        write(200,222) ncount*dt, sum(pop_wigner(1::2,ncount)), sum(pop_wigner(2::2,ncount))
     enddo

         
     write(*,*) "number of steps",nstep
     write(*,*) "ntraj",ntraj
     !call date_and_time(date,time,zone,timearray2)
     !write(*,*)"end ",timearray2
     !write(*,*)"diff ",timearray2-timearray
     !         write(*,*) "calculated energy", avg_energy
     !end if
         
222  format(200(e13.6,2x))
     !call MPI_Finalize(ierr)
     
     stop
   end program rpmd
!----------------------------------------------
subroutine diasym(eig,a,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','L',n,a,n,eig,work,l,inf)

end subroutine diasym
!-------------------------------------------------
