subroutine gethel(qs,hel)
  use global
  implicit none

  integer i,j,k,l,m,n

  real*8, intent(in)::  qs(ndof)
  real*8, intent(out):: hel(nstate,nstate)

  real*8, allocatable :: hmol(:,:)

  real*8 :: g,del
  
  allocate(hmol(molstate,molstate))
  !----------- definition of molecular hamiltonian--------------------
  hmol(1,2) = delta
  hmol(2,1) = hel(1,2)
  
  hmol(1,1) = 0.5*f0*qs(1)**2
  hmol(2,2) = 0.5*f0*(qs(1) - (shift))**2 - (bias)
  !==========================================================
  g = 0.005/27.2114

  !hel = 0.0d0

  ! diagonal elements
  do i=1,nstate
  do j=1,nstate
    m = int((i-1)/2) ! photon index 
    n = int((j-1)/2) ! photon index
    k = mod(i-1,2) ! 0 - D || 1 - A
    l = mod(j-1,2)
    hel(i,j) = hmol(k+1,l+1) * del(n,m) + (n + 0.5) * (omega_c) * del(i,j)  &
    + (g* sqrt(real(m+1)) * del(n,m+1) +  g* sqrt(real(m)) * del(n,m-1))*(1-del(k,l))
  end do
  end do

  return
end subroutine gethel
!==================================
Real*8 function del(i,j)

implicit none

Integer i,j

if(i.eq.j)then
 del = 1.0
else
 del = 0.0
end if

Return
End function del
!----------------------------