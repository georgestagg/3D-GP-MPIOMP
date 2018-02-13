subroutine initCond
    use workspace
    implicit none
    integer :: i, j ,k

    if(initialCondType .eq. 0) then
        GRID = 1.0d0
    else if (initialCondType .eq. 1) then
        call makeRandomPhase
    else if (initialCondType .eq. 2) then
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                    GRID(i,j,k) = 1.0d0*EXP(-EYE*atan2(GY(j)-NY/2*DSPACE,GX(i)-NX/2*DSPACE-2.0))&
                    *EXP(EYE*atan2(GY(j)-NY/2*DSPACE,GX(i)-NX/2*DSPACE+2.0))
            end do
        end do
    end do
    !$OMP end parallel do
    else
        GRID=1.0d0
    end if
end subroutine

subroutine loadPreviousState
    use params
    use output
    implicit none
    call read_wf_file(ICRfilename)
end subroutine

subroutine eulerStepOmega
    use params
    implicit none
    OMEGA = OMEGA + DOMEGADT*dble(DT)
    if(OMEGA < 0.0d0) then
        OMEGA = 0.0d0
    end if
end subroutine

subroutine makeRandomPhase
    use workspace
    implicit none
    double precision :: rpKC, rpAMP, phi
    integer :: i, j, k, n
    integer, dimension(:), allocatable :: seed

    call RANDOM_SEED(size = n)
    allocate(seed(n))
    seed = RANK
    call RANDOM_SEED(PUT = seed)

    call calc_rpKC_rpAMP(rpKC, rpAMP)
    if(RANK .eq. 0) then
        write (6, *) 'Imposing the highly non-equilibrium state...'
        write (6, *) 'K-space amplitude = ', rpAMP
        write (6, *) 'Maximum wavenumber = ', sqrt(rpKC)*DKSPACE
    end if

   !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                if (KX(i)**2+KY(j)**2+KZ(k)**2 <= rpKC*DKSPACE**2) then
                    call random_number(phi)
                    GRID(i,j,k) = rpAMP*exp(2.0*PI*EYE*phi)
                else
                    GRID(i,j,k) = 0.0d0
                end if
            end do
        end do
    end do
    !$OMP end parallel do

    call fftw_mpi_execute_dft(fftw_backward_plan, GRID, GRID)
    GRID=GRID/sqrt(dble(NX*NY*NZ))

end subroutine

subroutine calc_rpKC_rpAMP(rpKC, rpAMP)
    use params
    implicit none
    double precision, intent(out) :: rpKC, rpAMP
    double precision :: ev
    ev = ((0.5*NX/PI)**2.0)*ENERV
    rpAMP = sqrt(((0.11095279993106999d0*NV**2.5d0)*(NX*NY*NZ)) / (ev**1.5d0))
    rpKC = (1.666666666666666d0*ev)/NV
end subroutine