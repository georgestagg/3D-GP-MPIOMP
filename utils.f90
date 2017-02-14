subroutine setupGXYZ
    use params
    implicit none
    integer :: i, j ,k
    do i = PSX,PEX
        GX(i) = (i-1)*DSPACE
    end do
    
    do j = PSY,PEY
        GY(j) = (j-1)*DSPACE
    end do
    
    do k = PSZ,PEZ
        GZ(k) = (k-1)*DSPACE
    end do
end subroutine

subroutine initCond
    use params
    implicit none

    if(initialCondType .eq. 0) then
        GRID = 1.00d0
        call calc_POT
    end if

    if(initialCondType .eq. 2) then
        call loadPreviousState
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
