subroutine getdhel(qs,dhel)
  use global
  use pimc
  implicit none

  integer i,j,k,l,m,n

  real*8, intent(in)::  qs(ndof)

  real*8, intent(out)::  dhel(nstate,nstate,ndof)

  real*8 :: del2

  real*8, allocatable :: dhel_mol(:,:,:)
  allocate(dhel_mol(molstate,molstate,ndof))

  !---------- definiton of molecular gradient--------

  dhel_mol(1,1,1) = f0*qs(1)
  dhel_mol(2,2,1) = f0*(qs(1) - (shift))

  dhel_mol(1,2,1) = 0.0d0
  dhel_mol(2,1,1) = 0.0d0
  !-------------------------------------------------
  !dhel=0.
  
  do i=1,nstate
  do j=1,nstate
   
    m = int((i-1)/2) ! photon index 
    n = int((j-1)/2) ! photon index
    k = mod(i-1,2) ! 0 - D || 1 - A
    l = mod(j-1,2)

    dhel(i,j,1) = dhel_mol(k+1,l+1,1)*del2(n,m)
  

  enddo
  enddo
  
  !===========================================
  return
end subroutine getdhel
!==================================
Real*8 function del2(i,j)

implicit none

Integer i,j

if(i.eq.j)then
 del2 = 1.0
else
 del2 = 0.0
end if

Return
End function del2
!----------------------------
