
      subroutine evolve_freerp_exact(xe,ve)
         use global
         use pimc
         implicit none


         real*8, intent(inout)::xe(ndof,nb), ve(ndof,nb)

         real*8  p(ndof,nb)

         p = ve * mnuc

         call freerp (ndof,p,xe)

         ve = p / mnuc

      end subroutine evolve_freerp_exact

