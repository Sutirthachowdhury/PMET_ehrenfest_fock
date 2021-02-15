program rpmd
     use global
     use pimc,only:nbox,population,mnuc 
     !use globran
     !use mpi_lib

     implicit none
         
     !character(len=8)::date
     !character(len=10)::time
     !character(len=5)::zone
     !integer :: timearray(8)
     !integer :: timearray2(8)


     integer :: npre_eq,redmatsize,myid,nproc,proctraj,ierr
     integer :: istep,iup,idn     
     integer,allocatable :: seed(:)
     integer seed_dimension,time
   
     integer i,j,k,l,m,n,ncount
     real*8, allocatable:: xe(:,:),ve(:,:),fnon(:,:)
     real*8, allocatable:: xe_mc(:,:),x0_mc(:,:),p0_mc(:,:)
     real*8, allocatable:: x0(:,:),p0(:,:)
     real*8, allocatable:: corr_x(:),corr_pop(:),pop_coherent(:,:),pop_estimator_B(:),pop_estimator_A(:)
     real*8, allocatable:: xe_avg0(:)
     real*8, allocatable:: xe_avg(:,:)
     real*8, allocatable:: pop_wigner(:,:)
!     real*8, allocatable:: inv1_mc(:,:)
     real*8 :: dp,presum1,presum2,pref_mc,partition,pop_init,newA,boxlen,xmin,xmax,eterm
     real*8 :: gran,rand
 !    real*8 :: inv1_mc(nstate,nstate)
     real*8 :: corr_pop0,corr_x0,avg_energy,pop_init_nu,pop_init_de,B_init,population_initial
     real*8 :: hel(nstate,nstate),hel0
     real*8 :: coeff1(nstate,nstate)
     real*8 :: inv1(nstate,nstate)
     real*8 :: newpop_1(nstate),newpop_2(nstate)
     complex*16, allocatable:: redmat_traj(:,:,:)
     complex*16, allocatable :: redmat_traj_2(:,:,:)
!     complex*16, allocatable:: redmat(:,:,:)
     integer, parameter:: filelength = 1000

     ! IO for MC input parameters
     call parameter_input
     
     ! set up all the parameter for spin-boson model
     ! call setup

            
     ! Initialize the RNG seed
     call random_seed(size=seed_dimension)
     
     allocate(seed(seed_dimension))
     !write(*,*) "seed_dimension=",seed_dimension
     do i=1,seed_dimension
        seed(i)=time()+3*i-1
     enddo
     call random_seed(put=seed)




         allocate(xe_mc(ndof,nb),x0_mc(nstate,nb),p0_mc(nstate,nb))

         allocate(xe(ndof,nb),ve(ndof,nb),fnon(ndof,nb)) 
         allocate(x0(nstate,nb),p0(nstate,nb))

         allocate(corr_x(nstep),corr_pop(nstep),pop_estimator_B(nstep),pop_estimator_A(nstep))
         allocate(pop_coherent(nstate,nstep))
         allocate(xe_avg0(ndof))

         allocate(redmat_traj(nstate,nstate,nstep))
         allocate(redmat_traj_2(nstate,nstate,nstep))
       !  allocate(redmat(nstate,nstate,nstep))
         allocate(xe_avg(ndof,nstep))
         allocate(pop_wigner(nstate,nstep))
     !    allocate(inv1_mc(nstate,nstate))

        
         corr_x=0.0d0
         corr_pop=0.0d0

         pop_coherent = 0.0d0
         redmat_traj=cmplx(0.0d0,0.0d0)
         redmat_traj_2=cmplx(0.0d0,0.0d0)

        ! initialize all varialbes to 0 /or some value

         call init_rp(xe_mc)

         !         call plot_potential
 
         avg_energy = 0.
         partition = 0.
         corr_pop0 = 0.
         corr_x0 = 0.
         population = 0.

         pop_wigner = 0.
         
         do j=1, ntraj
            !call init_rp(xe_mc)
            call monte_sample_photo(xe_mc)
            call init_mapping_sampling(x0_mc,p0_mc)
                   
            xe=xe_mc
            x0=x0_mc
            p0=p0_mc
            
         
            write(10000 + int((j-1)/filelength)+1,222) real(j),((xe(l,k),k=1,nb),l=1,ndof),(x0(1,k),k=1,nb),&
                 (x0(2,k),k=1,nb),(p0(1,k),k=1,nb),(p0(2,k),k=1,nb)


            ! write(100,222) real(j),(xe(1,k),k=1,nb)
          ! write(101,222) real(j),(xe(2,k),k=1,nb)
          ! write(102,222) real(j),(xe(3,k),k=1,nb)
          ! write(104,222) real(j),(xe(60,k),k=1,nb)

           !-------------velocity distribution--------------------------------     
           
!            do k=1,ndof
!               dp=sqrt(omega(k)/(2.*tanh(0.5*beta*omega(k))))
!               dp = dsqrt(1.0/betap) !use reduced units   
!               do l=1,nb
!                  ve (k,l) = gran()*dp
!               end do
!            end do
            
           ! write(*,*) 've=',ve

!            ve = ve/mnuc
       !-----------------------------------------------------------------     
            ! operator x's value at t = 0.
!            xe_avg0 = 0.0d0
!            do i=1,nb
!               xe_avg0(:) = xe_avg0(:) + xe(:,i)
!            enddo
!            xe_avg0 = xe_avg0/real(nb)
 

            !============== OPERATOR B VALUE ===================
  !          B_init = 0.5*sum(x0(1,:)**2+p0(1,:)**2)/real(nb)-0.5

            !================================================================================================================            
            ! initial population value at time t = 0

 !           population_initial = 0.
            
 !           do i = 1,nb
 !              population_initial = population_initial + (0.5*(x0_mc(1,i)**2+p0_mc(1,i)**2-1.))
 !           enddo
          
 !           population_initial = population_initial/real(nb)

            !---------------- main dynamics propagation loop----------------------------

 !           call mapping_pot(xe,x0,p0,fnon)
             
 !           do istep=1,nstep
          
 !              call run_traj(xe,ve,x0,p0,fnon)
               
!               write(*,*) "fnon=",fnon(1,1),fnon(2,1),fnon(60,1)

             !========================= correlation picture ========================================================================== 
               ! for A-B case 
              ! corr_pop(istep) = corr_pop(istep) + pop_init*(0.5*sum(x0(1,:)**2+p0(1,:)**2)/real(nb)-0.5)*(pref_mc/dabs(pref_mc)) 
              
               !for B-B case
               !corr_pop(istep) = corr_pop(istep) + (B_init)*(0.5*sum(x0(1,:)**2+p0(1,:)**2)/real(nb)-0.5)*(pref_mc/dabs(pref_mc))
              
               !for A-A case
               !call updated_A(xe,ve,x0,p0,newA)
               !corr_pop(istep) = corr_pop(istep) +  (pop_init*newA)*(pref_mc/dabs(pref_mc))

               !for population-estimator
               !pop_estimator_B(istep) = pop_estimator_B(istep) + (0.5*sum(x0(1,:)**2+p0(1,:)**2)/real(nb)-0.5)*(pref_mc/dabs(pref_mc))
               !for population-estimator_A
               !pop_estimator_A(istep) = pop_estimator_A(istep) + (newA*(pref_mc/dabs(pref_mc)))
               
               !call population_nandini(x0,p0,newpop)
              !=============================================== end of correlation =======================================================

  !          do k=1,nstate
  !             do i=1,nb
  !                pop_wigner(k,istep) =  pop_wigner(k,istep) + (0.5*(x0(k,i)**2+p0(k,i)**2-1.0))/real(nb)
  !             enddo
  !          enddo

            
              
  !          do m=1,nstate
  !             do n=1,nstate  
  !                do i=1,nb
  !                   redmat_traj(m,n,istep) = redmat_traj(m,n,istep) & 
  !                        + (0.5*(x0_mc(1,i)**2+p0_mc(1,i)**2-1.0))*(0.5*(x0(m,i)+eye*p0(m,i))*(x0(n,i)-eye*p0(n,i)) & 
  !                        -0.5*del(m,n))/real(nb) 
  !                enddo
  !             enddo
  !          end do
          

  !            do m=1,nstate
  !             do n=1,nstate
  !                do i=1,nb
  !                   redmat_traj_2(m,n,istep) = redmat_traj_2(m,n,istep) &
  !                        + (0.5*(x0(m,i)+eye*p0(m,i))*(x0(n,i)-eye*p0(n,i))-0.5*del(m,n))/real(nb)
  !                        
  !                enddo
  !             enddo
  !          end do

               
           !end of time step loop
  !       enddo
         ! end of trajectory loop
         end do
        

         
222      format(1000(e13.6,2x))
             !call MPI_Finalize(ierr)
         
         stop
  end program
