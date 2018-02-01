program gp
    use workspace
    implicit none
    call init_params
    call init_parallel(RHSType)
    call parallel_barrier

    if(RANK .eq. 0) then
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "This is 3D-GP-MPIOMP - Written by GWS "
        write(6,'(a)') "Web: http://mas-gitlab.ncl.ac.uk/ngs54/            "
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Main parameters:"
        write(6,'(a,i5,a,i5,a,i5,a)') "NX: ", NX,     ", NY: ", NY    , ", NZ: ",  NZ, "."
        write(6,'(a,e10.3,a,e10.3,a)') "DSPACE: ", DSPACE, ", DTSIZE: ", DTSIZE, "."
        write(6,'(a,i8,a,i8,a)') "DUMPWF: ", DUMPWF, ", DUMPUTIL: ", DUMPUTIL, "."
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Initialising system..."
    end if
    call parallel_barrier

    call init_workspace
    call initCond
    call calc_POT
    call parallel_barrier
    call final_checks
    
    if(RANK .eq. 0) then
        write(6,'(a)') "Finished initialising!"
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Starting simulation..."
    end if
    call parallel_barrier

    DT = -EYE*DTSIZE
    call simulate(ISTEPS,0)
    DT = DTSIZE
    call simulate(NSTEPS,1)
    
    call parallel_barrier
    if(RANK .eq. 0) then
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Simulation all done!"
        write(6,'(a)') "---------------------------------------------------"
    end if
    call finalize_parallel
end PROGRAM gp

subroutine final_checks
    use output
    use rhs_RK4
    use rhs_FFTW
    implicit none
    if (METHOD==1) then
        call run_checks_FFTW
    end if
end subroutine

subroutine simulate(steps,rt)
    use output
    use rhs_RK4
    use rhs_FFTW
    implicit none
    integer :: steps,rt,i
    double precision :: perc

    do i = INITSSTEP, steps-1
        !Housekeeping
        if (modulo(i,DUMPUTIL) == 0) then
            if(RANK == 0) then
                if(steps == 0) then
                    perc = 100
                else
                    perc = dble(i)/steps*100.0d0
                end if
                if (rt == 1) then
                    write (6,fmt="(a,f6.2,a)") "Simulating: ",perc,"%"
                else
                    write (6,fmt="(a,f6.2,a)") "Ground State: ",perc,"%"
                end if
            end if
        end if

        if (modulo(i,DUMPWF) == 0) then
            call dump_wavefunction(i,rt)
        end if

        if(recalculatePot .and. rt == 1) then
            call calc_POT
        end if

        !Time stepping routines
        if(METHOD==0) then
            call RK4_step(rt)
        else if (METHOD==1) then
            call FFTW_step(rt)
        end if
        call eulerStepOmega
        TIME = TIME + dble(DT)
    end do
end subroutine