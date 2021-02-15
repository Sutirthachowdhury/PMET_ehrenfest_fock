subroutine freerp (nf,p,q)
      use global
      implicit real*8  (a-h,o-z)
!      integer*8 k,j,init
!     ------------------------------------------------------------------
!     Free harmonic ring-polymer evolution through a time interval dt.
!     ------------------------------------------------------------------
!
      dimension p(nf,nb),q(nf,nb)
      parameter (nbmax = 1024)
      dimension poly(4,nbmax)
      data init /0/
      save init,poly
!
      if (init .eq. 0) then
         if (nb .gt. nbmax) stop 'freerp 1'
         call ring (nf,poly)
         init = 1
      endif
      if (nb .eq. 1) then
         do j = 1,nf
            q(j,1) = q(j,1)+p(j,1)*poly(3,1)
         enddo
      else
         call realft (p,nf,nb,+1)
         call realft (q,nf,nb,+1)
         do k = 1,nb
            do j = 1,nf
               pjknew = p(j,k)*poly(1,k)+q(j,k)*poly(2,k)
               q(j,k) = p(j,k)*poly(3,k)+q(j,k)*poly(4,k)
               p(j,k) = pjknew
            enddo
         enddo

         call realft (p,nf,nb,-1)
         call realft (q,nf,nb,-1)
      endif

      return
      end subroutine

