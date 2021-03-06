module parallel
  use parallel_3DWithGhost
  use parallel_FFTW
  integer :: NPROCS, RANK, METHOD, MPI_WORLD, IERR, MPI_COMM
contains

  character(4096) function parallel_env_name()
    implicit none
    if (METHOD == 0) then
      parallel_env_name = "MPI/OMP with 3D block distribution."
    else if (METHOD == 1) then
      parallel_env_name = "MPI/OMP with 1D block distribution provided by FFTW3."
    end if
  end function

  subroutine init_parallel(RHSType)
    implicit none
    integer, intent(in) :: RHSType
    include 'mpif.h'
    METHOD = 0
    if (RHSType == 2) then
      METHOD = 1
    end if

    if (METHOD == 0) then
      call init_parallel_3DWithGhost
    else if (METHOD == 1) then
      call init_parallel_FFTW
    end if

    call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERR)

    if (RANK .eq. 0) then
      write (6, '(a)') "---------------------------------------------------"
      write (6, '(a,a)') "Initialised parallel environment: ", TRIM(parallel_env_name())
      write (6, '(a,i4,a)') "We're running ", NPROCS, " processes."
    end if
    call flush
    call parallel_barrier
    call run_omp_checks
    call setup_parallel_topology

    if (METHOD == 0) then
      MPI_COMM = MPI_COMM_GRID
    else if (METHOD == 1) then
      MPI_COMM = MPI_COMM_FFTW
    end if

  end subroutine

  subroutine run_omp_checks
    implicit none
    integer :: tid, nth
    !$OMP PARALLEL PRIVATE(tid,nth)
    nth = omp_get_num_threads()
    tid = omp_get_thread_num()
    if (tid .eq. 0) then
      write (6, '(a,i4,a,i4,a)') 'I am rank: ', RANK, ', I have ', nth, ' threads.'
    end if
    !$OMP END PARALLEL
  end subroutine

  subroutine setup_parallel_topology
    implicit none
    if (METHOD == 0) then
      call setup_parallel_topology_3DWithGhost(NPROCS, RANK)
    else if (METHOD == 1) then
      call setup_parallel_topology_FFTW
    end if
  end subroutine

  subroutine finalize_parallel
    implicit none
    if (METHOD == 0) then
      call finalize_parallel_3DWithGhost
    else if (METHOD == 1) then
      call finalize_parallel_FFTW
    end if
  end subroutine

  subroutine parallel_barrier
    implicit none
    call MPI_BARRIER(MPI_COMM, IERR_3DWG)
  end subroutine

end module
