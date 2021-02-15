module global
  implicit none

  integer, parameter:: nosc=1
  integer, parameter :: molstate = 2
  integer, parameter :: nfock = 5
  integer, parameter:: nstate=2*nfock
 
  integer,parameter:: ndof=nosc
  integer, parameter :: nbath = 100

  complex*16, allocatable:: redmat_sum(:,:,:)
 
  real*8, allocatable:: FREQ(:,:),TRANS(:,:,:)
  real*8, allocatable:: EVECTOR(:,:),eigen(:) !hel1(:,:) 

  real*8 omega(ndof),cq(nstate,ndof) ! harmonic oscillator frequency,coupling

  ! i.e. number of beads
  integer nb

  real*8 masse,beta,betap,dt,mom0,pos0,f0,shift,bias
  ! masse is the nuc mass ##

  real*8 omega_b(nbath), ctb(nbath)

  integer nstep,ntraj

  complex*16,parameter :: eye=(0.,1.)

  real*8 alpha,kappa,delta

  real*8 :: omega_c

  real*8  krondel(nstate,nstate) ! delta symbol 

  real*8 xi,wc
  integer :: ourseed

  real*8, parameter:: pi=3.1415926  
  real*8, parameter:: hbar=1.
end module global
