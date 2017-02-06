program gp
    use params
    use parallel
    use output
    implicit none
    call init_params
    call init_parallel
    call setup_parallel_topology
    call run_parallel_checks
    call calc_local_idx(NX,NY,NZ,PSX,PEX,PSY,PEY,PSZ,PEZ)

    if(RANK .eq. 0) then
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "This is 3D-GP-MPIMP - Written by GWS "
        write(6,'(a)') "Web: http://mas-gitlab.ncl.ac.uk/ngs54/            "
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Main parameters:"
        write(6,'(a,i5,a,i5,a,i5,a)') "NX: ", NX,     ", NY: ", NY    , ", NZ: ",  NZ, "."
        write(6,'(a,e10.3,a,e10.3,a)') "DSPACE: ", DSPACE, ", DTSIZE: ", DTSIZE, "."
        write(6,'(a,i8,a,i8,a)') "DUMPWF: ", DUMPWF, ", DUMPUTIL: ", DUMPUTIL, "."
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Initialising system geometry..."
    end if
    call MPI_BARRIER(COMM_GRID, IERR)
    call initialise
    call MPI_BARRIER(COMM_GRID, IERR)
    if(RANK .eq. 0) then
        write(6,'(a)') "Finished initialising!"
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Starting simulation..."
    end if
    DT = -EYE*DTSIZE
    call simulate(ISTEPS,0)
    DT = DTSIZE
    call simulate(NSTEPS,1)
    
    call MPI_BARRIER(COMM_GRID, IERR)
    if(RANK .eq. 0) then
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "All done!"
        write(6,'(a)') "---------------------------------------------------"
    end if
    call finalize_parallel()
end PROGRAM gp

subroutine initialise
    use params
    use parallel
    implicit none
    call init_arrays
    call setupGXYZ
    GRID = 0.0d0
    TIME = 0.0d0
    call initCond
end subroutine

subroutine simulate(steps,rt)
    use params
    use rhs
    use output
    use parallel
    implicit none
    integer :: steps,rt,i
    double precision :: perc

    !Timestep
    do i = INITSSTEP, steps-1
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
            if (rt == 1) then
                call dump_wavefunction(i)
            end if
        end if

        call RK4_step(rt)
        OMEGA = OMEGA + DOMEGADT*dble(DT)
        if(OMEGA < 0.0d0) then
            OMEGA = 0.0d0
        end if
        TIME = TIME + dble(DT)
        !if(potRep == 1 .and. rt == 1) then
        !    call calc_POT
        !end if
    end do
end subroutine