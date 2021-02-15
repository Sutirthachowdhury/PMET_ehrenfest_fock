subroutine Morse(x, qq, mm,uu, vv, mval)
  implicit none
  real*8, intent(in)  :: x,qq,mm,uu,vv
  real*8, intent(out) :: mval
  
  mval = (qq*(exp(-2.0*mm*(x-uu))-2*exp(-mm*(x-uu))) + vv)
  
  return
end subroutine Morse
!----------------------------------------------------------------
