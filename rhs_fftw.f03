module rhs_FFTW
  use workspace
contains
  subroutine run_checks_FFTW
    implicit none
    if (BCX == 0 .or. BCY == 0 .or. BCZ == 0) then
      if (RANK .eq. 0) then
        write (6, '(a)') "Warning: Reflective boundary conditions are not supported using the split-step fourier method."
        write (6, '(a)') "Forcing periodic boundary conditions..."
        BCX = 1
        BCY = 1
        BCZ = 1
      end if
    end if
    if (.not. (NX == NY .and. NY == NZ)) then
      if (RANK .eq. 0) then
        write (6, '(a)') "Error: Only NX=NY=NZ is currently implemented with the split-step fourier method."
        call finalize_parallel
        call exit(1)
      end if
    end if
    if (GAMMAC > 0.0d0) then
      if (RANK .eq. 0) then
        write (6, '(a)') "Warning: The dissipative GPE is not supported using the split-step fourier method."
        write (6, '(a)') "Forcing GAMMAC=0.0d0"
        GAMMAC = 0.0d0
      end if
    end if
    if (OMEGA > 0.0d0) then
      if (RANK .eq. 0) then
        write (6, '(a)') "Warning: A rotating frame is not supported using the split-step fourier method."
        write (6, '(a)') "Forcing OMEGA=0.0d0"
        OMEGA = 0.0d0
      end if
    end if
    if (FLUIDS > 1) then
      if (RANK .eq. 0) then
        write (6, '(a)') "Warning: Multi-component superfluid not supported when using the split-step fourier method."
        write (6, '(a)') "I'll just use the first one..."
      end if
    end if
  end subroutine

  subroutine FFTW_step
    implicit none
    integer:: i, j, k

    ! Half Step 1
    TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                  *conjg(WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez))
    call fftw_mpi_execute_dft(fftw_forward_plan_dens, TMPWS(1)%FLUID(1)%GRID, TMPWS(1)%FLUID(1)%GRID)
    TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                  /(NX*NY*NZ)*WS%FLUID(1)%DDI_K(sx:ex, sy:ey, sz:ez)
    call fftw_mpi_execute_dft(fftw_backward_plan_dens, TMPWS(1)%FLUID(1)%GRID, TMPWS(1)%FLUID(1)%GRID)

    WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                            *exp(-0.5d0*DT*EYE*(POT(sx:ex, sy:ey, sz:ez) &
                                                                + EDD_T1*WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                                *conjg(WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez)) &
                                                                + EDD_T2*TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez)))

    !Full step
    call fftw_mpi_execute_dft(fftw_forward_plan, WS%FLUID(1)%GRID, WS%FLUID(1)%GRID)
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz, ez
      do j = sy, ey
        do i = sx, ex
          WS%FLUID(1)%GRID(i, j, k) = WS%FLUID(1)%GRID(i, j, k)/(NX*NY*NZ) &
                                      *exp(-DT*EYE*0.5d0*(KX(i)**2.0d0 + KZ(k)**2.0d0 + KY(j)**2.0d0))
        end do
      end do
    end do
    !$OMP end parallel do
    call fftw_mpi_execute_dft(fftw_backward_plan, WS%FLUID(1)%GRID, WS%FLUID(1)%GRID)

    !Half Step 2
    TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                  *conjg(WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez))
    call fftw_mpi_execute_dft(fftw_forward_plan_dens, TMPWS(1)%FLUID(1)%GRID, TMPWS(1)%FLUID(1)%GRID)
    TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                  /(NX*NY*NZ)*WS%FLUID(1)%DDI_K(sx:ex, sy:ey, sz:ez)
    call fftw_mpi_execute_dft(fftw_backward_plan_dens, TMPWS(1)%FLUID(1)%GRID, TMPWS(1)%FLUID(1)%GRID)

    WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                            *exp(-0.5d0*DT*EYE*(POT(sx:ex, sy:ey, sz:ez) &
                                                                + EDD_T1*WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                                *conjg(WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez)) &
                                                                + EDD_T2*TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez)))
    if (RT == 0) then
      call FFTW_renormalise
    end if
  end subroutine

  subroutine FFTW_renormalise
    implicit none
    include 'mpif.h'
    double precision :: local_norm, total_norm, local_pot_int, total_pot_int

    TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) &
                                                  *conjg(WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez))

    local_norm = sum(TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez))*DSPACE*DSPACE*DSPACE
    call MPI_Allreduce(local_norm, total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, IERR)

    TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = 0.0d0
    where (POT < 1.0d0) TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = 1.0d0

    local_pot_int = sum(TMPWS(1)%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez))*DSPACE*DSPACE*DSPACE
    call MPI_Allreduce(local_pot_int, total_pot_int, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, IERR)

    WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez) = WS%FLUID(1)%GRID(sx:ex, sy:ey, sz:ez)/sqrt(total_norm)*sqrt(total_pot_int)

  end subroutine
end module
