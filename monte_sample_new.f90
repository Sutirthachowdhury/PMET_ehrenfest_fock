subroutine monte_sample_new(xe_mc)
  
  use pimc
  use global

  implicit none

  integer :: ibead,idof
  integer :: init_nmc
  integer :: i,j,iup,idn,k,init_mc

  real*8 :: iacc_x,iatt_x,rnum
  real*8 :: wtt,wt1,wt

  real*8, allocatable :: x1(:),x2(:),xnew(:,:)

  real*8, intent(inout) :: xe_mc(ndof,nb)
  
  !random number used in mc

  real*8 :: rand

  allocate(x1(ndof),x2(ndof),xnew(ndof,nb))
  
  !initialize all zeroes

  iatt_x = 0.d0
  iacc_x = 0.d0

  do init_mc=1,nmc
     !loop over all beads in each mc step                             
     do j = 1,ndof

        !randomely pick one particle
        call random_number(rand)
        idof = int(rand*ndof)+1
       ! write(*,*) "idof=", idof

        do i=1,nb

           !randomly pick a bead to move                                 
                                                
           call random_number(rand)
           ibead = int(rand*nb)+1
        !    write(*,*) "ibead=", ibead

           iup=ibead+1
           if(iup.gt.nb) iup=1
           idn=ibead-1
           if(idn.lt.1) idn=nb


           !move  nuclear coord
        
           !generating trial moves/weights
           
           ! do j=1,ndof
           call random_number(rand)
           rnum=rand
           xnew(idof,ibead) = xe_mc(idof,ibead)+(rnum-0.5d0)*step(1)
           x1(idof)=xe_mc(idof,idn);x2(idof)=xe_mc(idof,iup)
  !         write(*,*) "xnew=",xnew,"x1=",x1,"x2=",x2
           ! end do
           
           call init_nuclear_sampling(x1,xnew(idof,ibead),x2,wtt)
         !  write(*,*) 'wtt=',wtt   
           iatt_x=iatt_x+1.0
           
        
           !sampling call key
           
           !checking against original weight
           x1(idof)=xe_mc(idof,idn);x2(idof)=xe_mc(idof,iup)
           call init_nuclear_sampling(x1,xe_mc(idof,ibead),x2,wt1)
 !          write(*,*) 'wt1=',wt1          
           
           !accepting/rejecting step
           wt=wtt/wt1
!           write(*,*) 'wt=',wt
        
        call random_number(rand)
        rnum=rand
        

        if(rnum.lt.wt) then
           
           xe_mc(idof,ibead)=xnew(idof,ibead)
           iacc_x=iacc_x+1.d0
           
        endif

        
        !end loop over nb beads
     end do
     
  end do
  !     write(200,*) real(init_mc),xe_mc(1,1),xe_mc(2,1),xe_mc(3,1)
  

  !end loop over nmc
end do

!  print *, 'acceptance rate',iacc_x/iatt_x

end subroutine monte_sample_new
  !----------------------------------------------------------------------------------------------------
  subroutine init_nuclear_sampling(x1,x2,x3,wtt)


  use global
  use pimc
 ! use globran

  implicit none

  real*8, intent(in) :: x1(ndof),x2(ndof),x3(ndof)
  real*8, intent(out) :: wtt
  real*8 :: const
  integer :: k,l
 
  const = (0.5d0*mnuc/betap)

  wtt = 1.0d0

  wtt=wtt*dexp(-const*sum((x1(:)-x2(:))**2+(x2(:)-x3(:))**2))

  ! state independent part

  wtt = wtt*dexp(-betap*sum(0.5*mnuc*omega(:)**2*x2(:)**2))
!  write(*,*) 'wtt in sub=',wtt
  
end subroutine init_nuclear_sampling
!---------------------------------------------------------------------------------------------------------------
