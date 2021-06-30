program gp
  use workspace
  use init
  implicit none
  call init_params
  call init_parallel(RHSType)
  call parallel_barrier

  if (RANK .eq. 0) then
    write (6, '(a)') "---------------------------------------------------"
    write (6, '(a)') "This is 3D-GP-MPIOMP - Written by GWS "
    write (6, '(a)') "Web: http://mas-gitlab.ncl.ac.uk/ngs54/            "
    write (6, '(a)') "---------------------------------------------------"
    write (6, '(a)') "Main parameters:"
    write (6, '(a,i5,a,i5,a,i5,a)') "NX: ", NX, ", NY: ", NY, ", NZ: ", NZ, "."
    write (6, '(a,e10.3,a,e10.3,a)') "DSPACE: ", DSPACE, ", DTSIZE: ", DTSIZE, "."
    write (6, '(a,i8,a,i8,a)') "DUMPWF: ", DUMPWF, ", DUMPUTIL: ", DUMPUTIL, "."
    write (6, '(a,i3,a,1000f10.3)') "Simulating ", FLUIDS, " fluids with interaction matrix: ", GG
    write (6, '(a)') "---------------------------------------------------"
    write (6, '(a)') "Initialising system..."
  end if
  call parallel_barrier

  call init_workspace
  call initCond
  call calc_POT
  call parallel_barrier
  call final_checks

  if (RANK .eq. 0) then
    write (6, '(a)') "Finished initialising!"
    write (6, '(a)') "---------------------------------------------------"
    write (6, '(a)') "Starting simulation..."
  end if
  call parallel_barrier

  DT = -EYE*DTSIZE
  if (ISTEPS > 0 .and. TIME < 1e-14) then
    call simulate(ISTEPS, 0)
  end if
  DT = DTSIZE
  call simulate(NSTEPS, 1)

  call parallel_barrier
  if (RANK .eq. 0) then
    write (6, '(a)') "---------------------------------------------------"
    write (6, '(a)') "Simulation all done!"
    write (6, '(a)') "---------------------------------------------------"
  end if
  call finalize_parallel
end PROGRAM gp

subroutine final_checks
  use rhs_RK4
  use rhs_FFTW
  implicit none
  if (METHOD == 1) then
    call run_checks_FFTW
  end if
  if (METHOD == 0) then
    call run_checks_RK4
  end if
end subroutine

subroutine simulate(steps, set_rt)
  use io
  use rhs_RK4
  use rhs_FFTW
  implicit none
  integer :: steps, set_rt
  double precision :: perc

  RT = set_rt

  do cur_step = INITSTEP, steps
    !Housekeeping
    if (modulo(cur_step, DUMPUTIL) == 0) then
      if (RANK == 0) then
        if (steps == 0) then
          perc = 100
        else
          perc = dble(cur_step)/steps*100.0d0
        end if
        if (rt == 1) then
          write (6, fmt="(a,f6.2,a)") "Simulating: ", perc, "%"
        else
          write (6, fmt="(a,f6.2,a)") "Ground State: ", perc, "%"
        end if
      end if
    end if

    if (modulo(cur_step, DUMPWF) == 0) then
      call dump_wavefunction(cur_step)
    end if

    if (recalculatePot .and. rt == 1) then
      call calc_POT
    end if

    !Time stepping routines
    if (METHOD == 0) then
      call RK4_step
    else if (METHOD == 1) then
      call FFTW_step
    end if
    if (rt .eq. 1 .and. ABS(DOMEGADT) > 1e-14) then
      call eulerStepOmega
    end if
    TIME = TIME + dble(DT)
  end do
end subroutine
