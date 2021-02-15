module pimc

      implicit none

      !*** constants
!      real*8 :: pi

      !*** for histogram
      integer :: nbox,nboxstate
      real*8, allocatable :: delta_r0(:),pop(:,:),population(:)

      !*** for sampling_function
      ! no: of timeslices ie beads
     ! integer :: nb
      ! nuclear mass
      real*8 :: mnuc, hist_sum
      ! 1/kT
      !real*8 :: beta, betap

      !***for monte_sampling
      !# of monte carlo steps
     ! real*8 :: nmc
      integer :: nmc
      !* monte carlo step size
      real*8,allocatable :: step(:)
     
      !***for potential
      !====Ananth & Miller========
      !real*8 :: k1,k2
      !real*8 :: r1,r2,r12
      !real*8 :: eps1,eps2
      !real*8 :: gamma,c
      !============================
      
      !====Ananth 2013,JCP paper============
      !real*8 :: k1,c
      !real*8 :: a,delta
 
      !*** eta
      real*8 :: ef_eta

      !*** for potentials
      real*8, allocatable :: inv(:,:),coeff2(:,:)

      !**for average energy
      real*8 :: avg_energy

      !*** number of states
      !integer :: nstate,ndof

end module pimc
