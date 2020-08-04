module rhs_RK4
  use workspace
  use derivs
contains
  subroutine run_checks_RK4
    implicit none
    if (initialCondType == 2) then
      if (RANK .eq. 0) then
        write (6, '(a)') "Error: Highly non-equilibrium initialCondType is not supported with 3D block distribution."
        call finalize_parallel
        call exit(1)
      end if
    end if
  end subroutine

  subroutine RK4_step
    implicit none
    integer :: f
    call halo_swap_WS(WS)
    call RK4_calc_WS_RHS(WS, TMPWS(1))

    call RK4_sub_step(0.5d0, 6.0d0, WS)

    call halo_swap_WS(TMPWS(2))
    call RK4_calc_WS_RHS(TMPWS(2), TMPWS(1))

    call RK4_sub_step(0.5d0, 3.0d0, TMPWS(3))

    call halo_swap_WS(TMPWS(2))
    call RK4_calc_WS_RHS(TMPWS(2), TMPWS(1))

    call RK4_sub_step(1.0d0, 3.0d0, TMPWS(3))

    call halo_swap_WS(TMPWS(2))
    call RK4_calc_WS_RHS(TMPWS(2), TMPWS(1))

    call RK4_final_step

    !Renormalise WF after imaginary time decay
    if (RT .eq. 0 .and. (.not. NORENORM)) then
      do f = 1, FLUIDS
        call RK4_renormalise_fluid(WS%FLUID(f)%GRID, TMPWS(1)%FLUID(f)%GRID)
      end do
    end if
  end subroutine

  subroutine RK4_sub_step(d1, d2, IWS)
    implicit none
    integer :: f, i, j, k
    double precision :: d1, d2
    type(workspace_t), intent(in) :: IWS
    !$OMP parallel do private (f,i,j,k) collapse(4)
    do f = 1, FLUIDS
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            TMPWS(2)%FLUID(f)%GRID(i, j, k) = WS%FLUID(f)%GRID(i, j, k) + TMPWS(1)%FLUID(f)%GRID(i, j, k)*d1*DT
            TMPWS(3)%FLUID(f)%GRID(i, j, k) = IWS%FLUID(f)%GRID(i, j, k) + TMPWS(1)%FLUID(f)%GRID(i, j, k)*(DT/d2)
          end do
        end do
      end do
    end do
    !$OMP end parallel do
  end subroutine

  subroutine RK4_final_step
    implicit none
    integer :: f, i, j, k
    !$OMP parallel do private (f,i,j,k) collapse(4)
    do f = 1, FLUIDS
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            WS%FLUID(f)%GRID(i, j, k) = TMPWS(3)%FLUID(f)%GRID(i, j, k) + TMPWS(1)%FLUID(f)%GRID(i, j, k)*(DT/6.0d0)
          end do
        end do
      end do
    end do
    !$OMP end parallel do
  end subroutine

  subroutine RK4_calc_WS_RHS(ws_in, ws_out)
    implicit none
    type(workspace_t), intent(in) :: ws_in
    type(workspace_t), intent(out) :: ws_out
    integer :: f
    do f = 1, FLUIDS
      call RK4_calc_RHS_Fluid(f, ws_in, ws_out)
    end do
  end subroutine

  subroutine RK4_calc_RHS_Fluid(f, ws_in, ws_out)
    implicit none
    integer :: f, i, j, k, p
    type(workspace_t), intent(in) :: ws_in
    type(workspace_t), intent(out) :: ws_out

    if (RHSType .eq. 0) then
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            ws_out%FLUID(f)%GRID(i, j, k) = -0.5d0*laplacian(ws_in%FLUID(f), i, j, k) &
                                            + POT(i, j, k)*ws_in%FLUID(f)%GRID(i, j, k) &
                                            - ws_in%FLUID(f)%GRID(i, j, k)
            if (VELX > 0) then
              ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) + VELX*EYE*ddx(ws_in%FLUID(f), i, j, k)
            end if
            if (VELY > 0) then
              ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) + VELY*EYE*ddy(ws_in%FLUID(f), i, j, k)
            end if
            if (VELZ > 0) then
              ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) + VELZ*EYE*ddz(ws_in%FLUID(f), i, j, k)
            end if
            if (OMEGA > 0) then
              ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) + OMEGA*EYE*(GX(i)*ddy(ws_in%FLUID(f), i, j, k) &
                                                                                         - GY(j)*ddx(ws_in%FLUID(f), i, j, k))
            end if
          end do
        end do
      end do
      !$OMP end parallel do
      do p = 1, FLUIDS
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz + NGHOST, ez - NGHOST
          do j = sy + NGHOST, ey - NGHOST
            do i = sx + NGHOST, ex - NGHOST
              ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) + GG(f, p)*ws_in%FLUID(p)%GRID(i, j, k) &
                                              *CONJG(ws_in%FLUID(p)%GRID(i, j, k))*ws_in%FLUID(f)%GRID(i, j, k)
            end do
          end do
        end do
        !$OMP end parallel do
      end do
    end if
    if (RHSType .eq. 1) then
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            ws_out%FLUID(f)%GRID(i, j, k) = -0.5d0*laplacian(ws_in%FLUID(f), i, j, k) &
                                            + harm_osc_C*ws_in%FLUID(f)%GRID(i, j, k)*ws_in%FLUID(f)%GRID(i, j, k) &
                                            *CONJG(ws_in%FLUID(f)%GRID(i, j, k)) + POT(i, j, k)*ws_in%FLUID(f)%GRID(i, j, k) &
                                            - harm_osc_mu*ws_in%FLUID(f)%GRID(i, j, k) &
                                            + OMEGA*EYE*(GX(i)*ddy(ws_in%FLUID(f), i, j, k) - GY(j)*ddx(ws_in%FLUID(f), i, j, k))
          end do
        end do
      end do
      !$OMP end parallel do
    end if
    if (RHSType .eq. 3) then
      TMPWS(4)%FLUID(1)%field_number = f
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz, ez
        do j = sy, ey
          do i = sx, ex
            TMPWS(4)%FLUID(1)%GRID(i, j, k) = CONJG(exp(EYE*WS%FLUID(f)%QP_E(i, j, k, 1)))*ws_in%FLUID(f)%GRID(i, j, k)
          end do
        end do
      end do
      !$OMP end parallel do
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            ws_out%FLUID(f)%GRID(i, j, k) = -ws_in%FLUID(f)%GRID(i, j, k) - exp(EYE*WS%FLUID(f)%QP_E(i, j, k, 1)) &
                                            *d2dx2(TMPWS(4)%FLUID(1), i, j, k)
          end do
        end do
      end do
      !$OMP end parallel do
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz, ez
        do j = sy, ey
          do i = sx, ex
            TMPWS(4)%FLUID(1)%GRID(i, j, k) = CONJG(exp(EYE*WS%FLUID(f)%QP_E(i, j, k, 2)))*ws_in%FLUID(f)%GRID(i, j, k)
          end do
        end do
      end do
      !$OMP end parallel do
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) - exp(EYE*WS%FLUID(f)%QP_E(i, j, k, 2)) &
                                            *d2dy2(TMPWS(4)%FLUID(1), i, j, k)
          end do
        end do
      end do
      !$OMP end parallel do
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz, ez
        do j = sy, ey
          do i = sx, ex
            TMPWS(4)%FLUID(1)%GRID(i, j, k) = CONJG(exp(EYE*WS%FLUID(f)%QP_E(i, j, k, 3)))*ws_in%FLUID(f)%GRID(i, j, k)
          end do
        end do
      end do
      !$OMP end parallel do
      !$OMP parallel do private (i,j,k) collapse(3)
      do k = sz + NGHOST, ez - NGHOST
        do j = sy + NGHOST, ey - NGHOST
          do i = sx + NGHOST, ex - NGHOST
            ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) - exp(EYE*WS%FLUID(f)%QP_E(i, j, k, 3)) &
                                            *d2dz2(TMPWS(4)%FLUID(1), i, j, k)
          end do
        end do
      end do
      !$OMP end parallel do

      do p = 1, FLUIDS
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz + NGHOST, ez - NGHOST
          do j = sy + NGHOST, ey - NGHOST
            do i = sx + NGHOST, ex - NGHOST
              ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k) + GG(f, p)*ws_in%FLUID(p)%GRID(i, j, k) &
                                              *CONJG(ws_in%FLUID(p)%GRID(i, j, k))*ws_in%FLUID(f)%GRID(i, j, k)
            end do
          end do
        end do
        !$OMP end parallel do
      end do
    end if

    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz + NGHOST, ez - NGHOST
      do j = sy + NGHOST, ey - NGHOST
        do i = sx + NGHOST, ex - NGHOST
          ws_out%FLUID(f)%GRID(i, j, k) = ws_out%FLUID(f)%GRID(i, j, k)/(EYE - GAMMAC)
        end do
      end do
    end do
    !$OMP end parallel do
  end subroutine

  subroutine halo_swap_WS(ws_in)
    implicit none
    type(workspace_t) :: ws_in
    integer :: f, n
    do f = 1, FLUIDS
      do n = 1, NGHOST
        call halo_swap_complex(ws_in%FLUID(f)%GRID, n)
      end do
    end do
  end subroutine

  subroutine halo_swap_complex(kk, n)
    use parallel
    implicit none
    complex*16, dimension(sx:ex, sy:ey, sz:ez) :: kk
    integer :: n
    include 'mpif.h'
    !X dim - tag 0 and 1
    call MPI_sendrecv(kk(ex - NGHOST - n + 1, sy:ey, sz:ez), (ey - sy + 1)*(ez - sz + 1), MPI_DOUBLE_COMPLEX, cart_shift(0)%rd, 0, &
                      kk(sx + NGHOST - n, sy:ey, sz:ez), (ey - sy + 1)*(ez - sz + 1), &
                      MPI_DOUBLE_COMPLEX, cart_shift(0)%rs, 0, MPI_COMM_GRID, MPISTAT, IERR)
    if (IERR .ne. MPI_SUCCESS) then
      write (6, *) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
      write (6, *) "Something has gone wrong... Quitting."
      CALL EXIT(1)
    end if
    call MPI_sendrecv(kk(sx + NGHOST + n - 1, sy:ey, sz:ez), (ey - sy + 1)*(ez - sz + 1), MPI_DOUBLE_COMPLEX, cart_shift(0)%rs, 1, &
                      kk(ex - NGHOST + n, sy:ey, sz:ez), (ey - sy + 1)*(ez - sz + 1), &
                      MPI_DOUBLE_COMPLEX, cart_shift(0)%rd, 1, MPI_COMM_GRID, MPISTAT, IERR)
    if (IERR .ne. MPI_SUCCESS) then
      write (6, *) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
      write (6, *) "Something has gone wrong... Quitting."
      CALL EXIT(1)
    end if
    !Y dim - tag 2 and 3
    call MPI_sendrecv(kk(sx:ex, ey - NGHOST - n + 1, sz:ez), (ex - sx + 1)*(ez - sz + 1), MPI_DOUBLE_COMPLEX, cart_shift(1)%rd, 2, &
                      kk(sx:ex, sy + NGHOST - n, sz:ez), (ex - sx + 1)*(ez - sz + 1), &
                      MPI_DOUBLE_COMPLEX, cart_shift(1)%rs, 2, MPI_COMM_GRID, MPISTAT, IERR)
    if (IERR .ne. MPI_SUCCESS) then
      write (6, *) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
      write (6, *) "Something has gone wrong... Quitting."
      CALL EXIT(1)
    end if
    call MPI_sendrecv(kk(sx:ex, sy + NGHOST + n - 1, sz:ez), (ex - sx + 1)*(ez - sz + 1), MPI_DOUBLE_COMPLEX, cart_shift(1)%rs, 3, &
                      kk(sx:ex, ey - NGHOST + n, sz:ez), (ex - sx + 1)*(ez - sz + 1), &
                      MPI_DOUBLE_COMPLEX, cart_shift(1)%rd, 3, MPI_COMM_GRID, MPISTAT, IERR)
    if (IERR .ne. MPI_SUCCESS) then
      write (6, *) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
      write (6, *) "Something has gone wrong... Quitting."
      CALL EXIT(1)
    end if
    !Z dim - tag 4 and 5
    call MPI_sendrecv(kk(sx:ex, sy:ey, ez - NGHOST - n + 1), (ex - sx + 1)*(ey - sy + 1), MPI_DOUBLE_COMPLEX, cart_shift(2)%rd, 4, &
                      kk(sx:ex, sy:ey, sz + NGHOST - n), (ex - sx + 1)*(ey - sy + 1), &
                      MPI_DOUBLE_COMPLEX, cart_shift(2)%rs, 4, MPI_COMM_GRID, MPISTAT, IERR)
    if (IERR .ne. MPI_SUCCESS) then
      write (6, *) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
      write (6, *) "Something has gone wrong... Quitting."
      CALL EXIT(1)
    end if
    call MPI_sendrecv(kk(sx:ex, sy:ey, sz + NGHOST + n - 1), (ex - sx + 1)*(ey - sy + 1), MPI_DOUBLE_COMPLEX, cart_shift(2)%rs, 5, &
                      kk(sx:ex, sy:ey, ez - NGHOST + n), (ex - sx + 1)*(ey - sy + 1), &
                      MPI_DOUBLE_COMPLEX, cart_shift(2)%rd, 5, MPI_COMM_GRID, MPISTAT, IERR)
    if (IERR .ne. MPI_SUCCESS) then
      write (6, *) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
      write (6, *) "Something has gone wrong... Quitting."
      CALL EXIT(1)
    end if
  end subroutine

  subroutine RK4_renormalise_fluid(gt, tgt)
    implicit none
    complex*16, dimension(sx:ex, sy:ey, sz:ez) :: gt, tgt
    include 'mpif.h'
    double precision :: local_norm, total_norm, local_pot_int, total_pot_int

    tgt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST) = &
      gt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST) &
      *conjg(gt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST))

    local_norm = sum(tgt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST))*DSPACE*DSPACE*DSPACE
    call MPI_Allreduce(local_norm, total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, IERR)

    tgt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST) = 0.0d0
    where (POT < 1.0d0) tgt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST) = 1.0d0

    local_pot_int = sum(tgt(sx + NGHOST:ex - NGHOST, sy + NGHOST:ey - NGHOST, sz + NGHOST:ez - NGHOST))*DSPACE*DSPACE*DSPACE
    call MPI_Allreduce(local_pot_int, total_pot_int, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW, IERR)

    gt(sx:ex, sy:ey, sz:ez) = gt(sx:ex, sy:ey, sz:ez)/sqrt(total_norm)*sqrt(total_pot_int)
  end subroutine
end module
