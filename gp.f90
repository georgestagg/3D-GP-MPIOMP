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
        write(6,'(a)') "This is 3D-GP-MPIMP - Original code written by GWS "
        write(6,'(a)') "Web: http://mas-gitlab.ncl.ac.uk/ngs54/            "
        write(6,'(a)') "---------------------------------------------------"
        write(6,'(a)') "Main parameters:"
        write(6,'(a,i5,a,i5,a,i5,a)') "NX: ", NX,     ", NY: ", NY    , ", NZ: ",  NZ, "."
        write(6,'(a,e10.3,a,e10.3,a)') "DSPACE: ", DSPACE, ", DTSIZE: ", DTSIZE, "."
        write(6,'(a)') "---------------------------------------------------"
    end if
    call MPI_BARRIER(COMM_GRID, IERR)
    call initialise
    DT = -EYE*DTSIZE
    !call runit(ISTEPS,0,0)
    !DT = DTSIZE
    !call runit(NSTEPS,1,1)
    
    call make_file("./test.dat")
    call write_wf_file
    call close_file
end PROGRAM gp

subroutine initialise
    use params
    use parallel
    call init_arrays
    call setupGXYZ
    call setupMeshgrid
    call calc_POT
    GRID = 0.0d0
    TIME = 0.0d0
    call initCond
end subroutine