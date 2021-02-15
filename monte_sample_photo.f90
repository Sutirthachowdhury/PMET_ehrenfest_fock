subroutine monte_sample_photo(xe_mc)
  
  use pimc
  use global

  implicit none

  integer :: ibead,imove
  integer :: init_nmc
  integer :: i,j,iup,idn,k,init_mc

  real*8 :: iacc_x,iatt_x,rnum
  real*8 :: energy,energy_old,del_e,wt

  real*8, allocatable :: x1(:),x2(:),xnew(:)

  real*8, intent(inout) :: xe_mc(ndof,nb)
  
  !random number used in mc

  real*8 :: rand

  allocate(x1(ndof),x2(ndof),xnew(ndof))
  
  !initialize all zeroes

  iatt_x = 0.d0
  iacc_x = 0.d0

  do init_mc=1,nmc
     !loop over all beads in each mc step                             
     
     !do i=1,nb

        !randomly pick a bead to move                                 
                                                
        call random_number(rand)
        ibead = int(rand*nb)+1
        ! write(*,*) "ibead=", ibead

        iup=ibead+1
        if(iup.gt.nb) iup=1
        idn=ibead-1
        if(idn.lt.1) idn=nb


        !move  nuclear coord
        
        !generating trial moves/weights
        
        do j=1,ndof
           call random_number(rand)
           rnum=rand
           xnew(j) = xe_mc(j,ibead)+(rnum-0.5d0)*step(1)
           x1(j)=xe_mc(j,idn);x2(j)=xe_mc(j,iup)
  !         write(*,*) "xnew=",xnew,"x1=",x1,"x2=",x2
        end do

        call init_nuclear_sampling(x1,xnew,x2,energy)
!        write(*,*) 'wtt=',wtt   
        iatt_x=iatt_x+1.0
           
        
        !sampling call key
        !checking against original weight
        
        
        x1(:)=xe_mc(:,idn);x2(:)=xe_mc(:,iup)
        call init_nuclear_sampling(x1,xe_mc(:,ibead),x2,energy_old)


        ! calculating the energy difference
        del_e = (energy-energy_old)

        !accepting/rejecting step
        wt=dexp(-betap*del_e)
        
        
        call random_number(rand)
        rnum=rand
        

        if(rnum.lt.wt) then
           
           xe_mc(:,ibead)=xnew(:)
           iacc_x=iacc_x+1.d0
           
        endif

        
        !end loop over nb beads
     !end do


!     write(200,*) real(init_mc),xe_mc(1,1),xe_mc(2,1),xe_mc(3,1)


     !end loop over nmc
     end do

  !print *, 'acceptance rate',iacc_x/iatt_x

end subroutine monte_sample_photo
  !----------------------------------------------------------------------------------------------------
  subroutine init_nuclear_sampling(x1,x2,x3,energy)


  use global
  use pimc
 ! use globran

  implicit none

  real*8, intent(in) :: x1(ndof),x2(ndof),x3(ndof)
  real*8, intent(out) :: energy
  real*8 :: const
  integer :: k,l
 
  const = ((0.5d0*mnuc)/(betap)**2)


  energy = const*((x1(ndof)-x2(ndof))**2+(x2(ndof)-x3(ndof))**2)  &
       + (0.5*mnuc*(omega(1)**2)*(x2(ndof) - 7.0)**2)
  
  !energy = const*((x1(ndof)-x2(ndof))**2)  &
  !     + (0.5*mnuc*(omega(1)**2)*(x2(ndof) - 7.0)**2)

end subroutine init_nuclear_sampling
!---------------------------------------------------------------------------------------------------------------
