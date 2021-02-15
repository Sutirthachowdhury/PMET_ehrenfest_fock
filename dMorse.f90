subroutine dMorse(x, qq, mm,uu, mval)
  implicit none
  real*8, intent(in)  :: x,qq,mm,uu
  real*8, intent(out) :: mval

  mval = (qq*(exp(-2.0*mm*(x-uu))*(-2.0*mm) + 2*exp(-mm*(x-uu))*mm))
  
  return
end subroutine dmorse
