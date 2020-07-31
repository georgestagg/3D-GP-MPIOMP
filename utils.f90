subroutine eulerStepOmega
  use params
  implicit none
  OMEGA = OMEGA + DOMEGADT*dble(DT)
  if (OMEGA < 0.0d0) then
    OMEGA = 0.0d0
  end if
end subroutine
