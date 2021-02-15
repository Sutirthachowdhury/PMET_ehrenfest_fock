 subroutine ring (nf,poly)
   use global, only:  masse,beta,hbar,nb,dt
   use pimc
   implicit real*8  (a-h,o-z)
!   integer*8  k
!     ------------------------------------------------------------------
!     Monodromy matrix elements for free ring-polymer evolution.
!     ------------------------------------------------------------------
!
     ! common /system/ beta,hbar,em
      dimension poly(4,nb)
!
      poly(1,1) = 1.d0
      poly(2,1) = 0.d0
      poly(3,1) = dt/mnuc
      poly(4,1) = 1.d0
      if (nb .gt. 1) then
         betan = beta/nb
         twown = 2.d0/(betan*hbar)
         pibyn = dacos(-1.d0)/nb
         do k = 1,nb/2
            wk = twown*dsin(k*pibyn)
            wt = wk*dt
            wm = wk*mnuc
            cwt = dcos(wt)
            swt = dsin(wt)
            poly(1,k+1) = cwt
            poly(2,k+1) = -wm*swt
            poly(3,k+1) = swt/wm
            poly(4,k+1) = cwt
         enddo
         do k = 1,(nb-1)/2
            poly(1,nb-k+1) = poly(1,k+1)
            poly(2,nb-k+1) = poly(2,k+1)
            poly(3,nb-k+1) = poly(3,k+1)
            poly(4,nb-k+1) = poly(4,k+1)
         enddo
      endif
      return
      end subroutine


