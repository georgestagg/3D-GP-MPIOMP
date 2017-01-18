subroutine calc_POT
    use params
    implicit none
    if (enableTrap) then
        POT = 0.5d0*(TXSCALE*(MGX-TX)**2.0d0+TYSCALE*(MGY-TY)**2.0d0+TZSCALE*(MGZ-TZ)**2.0d0)
    end if
end subroutine