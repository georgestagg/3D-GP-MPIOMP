module parallel
    use parallel_3DWithGhost
    use parallel_FFTW
    integer :: NNODES,RANK,METHOD,MPI_WORLD,IERR
    contains

    character(4096) function parallel_env_name()
        implicit none
        if(METHOD==0) then
            parallel_env_name = "MPI/OMP with 3D block distribution."
        else if(METHOD==1) then
            parallel_env_name = "MPI/OMP with 1D block distribution provided by FFTW3."
        end if
    end function

    subroutine init_parallel(RHSType)
        implicit none
        integer :: MPI_WORLD
        integer,intent(in) :: RHSType
        include 'mpif.h'
        if(RHSType < 2) then
            METHOD = 0
        else if (RHSType .eq. 2) then
            METHOD = 1
        end if

        MPI_WORLD = MPI_COMM_WORLD

        if(METHOD==0) then
            call init_parallel_3DWithGhost
        else if(METHOD==1) then
            call init_parallel_FFTW
        end if

        call MPI_COMM_RANK(MPI_WORLD, RANK, IERR)
        call MPI_COMM_SIZE(MPI_WORLD, NNODES, IERR)

        if(RANK .eq. 0) then
            write(6,'(a)') "---------------------------------------------------"
            write(6,'(a,a)') "Initialised parallel environment: ", TRIM(parallel_env_name())
            write(6,'(a,i4,a)') "We're running on ",NNODES, " nodes."
        end if
        call flush
        call parallel_barrier
        call run_omp_checks
        call setup_parallel_topology
        
    end subroutine

    subroutine run_omp_checks
        implicit none
        integer :: tid, nth
        !$OMP PARALLEL PRIVATE(tid,nth)
            nth = omp_get_num_threads()
            tid = omp_get_thread_num()
            if(tid .eq. 0) then
                write(6, '(a,i4,a,i4,a)') 'I am rank: ', RANK,', I have ',nth,' threads.'
            end if
        !$OMP END PARALLEL
    end subroutine

    subroutine setup_parallel_topology
        implicit none
        if(METHOD==0) then
            call setup_parallel_topology_3DWithGhost(NNODES,RANK)
        else if(METHOD==1) then
            call setup_parallel_topology_FFTW
        end if
    end subroutine


    subroutine finalize_parallel
        implicit none
        if(METHOD==0) then
            call finalize_parallel_3DWithGhost
        else if(METHOD==1) then
            call finalize_parallel_FFTW
        end if
    end subroutine

    subroutine parallel_barrier
        implicit none
        if(METHOD==0) then
            call MPI_BARRIER(COMM_GRID, IERR_3DWG)
        else if(METHOD==1) then
            call MPI_BARRIER(MPI_WORLD, IERR)
        end if
    end subroutine
end module