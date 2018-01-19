subroutine initCond
    use workspace
    implicit none

    if(initialCondType .eq. 0) then
        GRID = 1.00d0
        call calc_POT
    end if
    
end subroutine

subroutine loadPreviousState
    use params
    use output
    implicit none
    call read_wf_file(ICRfilename)
    INITSSTEP = RESUMESTEP
    TIME = RESUMETIME
end subroutine

subroutine eulerStepOmega
    use params
    implicit none
    OMEGA = OMEGA + DOMEGADT*dble(DT)
    if(OMEGA < 0.0d0) then
        OMEGA = 0.0d0
    end if
end subroutine