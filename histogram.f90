subroutine histogram(xe_mc,term)

  use pimc,only : nbox,population
  use global,only : nb,nstate,ndof,betap

  real*8, intent(in) :: xe_mc(ndof,2*nb)
  real*8, intent(in) :: term

  integer :: ibead, ibox, k, i
  real*8 :: boxlen,xmin,xmax,hist_sum

  xmin = -10.0d0
  xmax = 10.0d0
  boxlen = (xmax-xmin)/nbox
 ! hist_sum = sum(population)

  !do ibead=1,nb
      ibox = int((xe_mc(ndof,nb)-xmin)/boxlen)+1
      !if x in range, returns from 1 to nbox
      if ((ibox.le.nbox).and.(ibox.gt.0)) then
         population(ibox) = population(ibox) + term
         !delta_r0(ibox) = delta_r0(ibox) + term
         !write(*,*) "term=",term
         !pop(:,ibox)=pop(:,ibox)+dpref(:)
      end if
      !hist_sum = sum(population)
      !if(hist_sum.gt.0.d0) population = population/hist_sum/boxlen
!      write(*,*) "delta",delta_r0
  !end do

end subroutine histogram
