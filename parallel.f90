module parallel
    integer :: NNODES,NODE_DIMS(3),NODE_COORDS(3)
    logical :: PERIODIC(3)
    integer :: COMM_WORLD,COMM_GRID,RANK,IERR
    integer, dimension(:), ALLOCATABLE :: MPISTAT
    contains
    subroutine init_parallel
        implicit none
        character(len=1024) :: check_NNODES_text
        include 'mpif.h'

        ALLOCATE(MPISTAT(MPI_STATUS_SIZE))
        COMM_WORLD = MPI_COMM_WORLD
        call MPI_INIT(IERR)
        call MPI_COMM_RANK(COMM_WORLD, RANK, IERR)
        call MPI_COMM_SIZE(COMM_WORLD, NNODES, IERR)

        if(RANK .eq. 0) then
            write(6,'(a)') "---------------------------------------------------"
            write(6,'(a)') "Init Parallel:"
            write(6,'(a,i4,a)') "We're running on ",NNODES, ' nodes.'
        end if

        call MPI_BARRIER(COMM_WORLD, IERR)
    end subroutine

    subroutine setup_parallel_topology
        implicit none
        !These are not GPE periodicity, it is the periodicity of the MPI topology
        PERIODIC(1) = .true.
        PERIODIC(2) = .true.
        PERIODIC(3) = .true.
        NODE_DIMS = 0
        call MPI_DIMS_CREATE(NNODES, 3, NODE_DIMS,IERR)
        call MPI_CART_CREATE(COMM_WORLD, 3, NODE_DIMS, PERIODIC, .true., COMM_GRID, IERR)
        if(RANK .eq. 0) then
            write(6,'(a,i3,a,i3,a,i3,a)') "Topology is: (", NODE_DIMS(1),",", &
                NODE_DIMS(2),",", NODE_DIMS(3),")."
        end if
        call MPI_COMM_RANK(COMM_GRID, RANK, IERR)
        call MPI_BARRIER(COMM_GRID, IERR)
    end subroutine

    subroutine run_parallel_checks
        use omp_lib
        implicit none
        integer :: tid, nth
        !$OMP PARALLEL PRIVATE(tid,nth)
            nth = omp_get_num_threads()
            tid = omp_get_thread_num()
            if(tid .eq. 0) then
                write(6, '(a,i4,a,i3,a)') 'I am rank: ',rank,', I have ',nth,' threads.'
            end if
        !$OMP END PARALLEL
        call MPI_BARRIER(COMM_GRID, IERR)
    end subroutine

    subroutine finalize_parallel
        implicit none
        call MPI_FINALIZE(IERR)
    end subroutine

    !get local grid sizes
    subroutine calc_local_idx(NX,NY,NZ,PSX,PEX,PSY,PEY,PSZ,PEZ)
        implicit none
        integer :: pnx,pny,pnz
        integer,intent(in)  :: NX,NY,NZ
        integer,intent(out) :: PSX,PEX,PSY,PEY,PSZ,PEZ
        pnx = NX/NODE_DIMS(1)
        pny = NY/NODE_DIMS(2)
        pnz = NY/NODE_DIMS(3)
        !if(RANK .eq. 0) then
        !    write(6,'(a,i5,a,i5,a,i5,a)') "PNX: ", pnx,     ", PNY: ", pny    , ", PNZ: ",  pnz, "."
        !end if
        
        call MPI_CART_COORDS(COMM_GRID,RANK,3,NODE_COORDS,IERR)

        !Get local grid (plus ghost)
        PSX = NODE_COORDS(1)*pnx
        PEX = NODE_COORDS(1)*pnx + pnx + 1
        PSY = NODE_COORDS(2)*pny
        PEY = NODE_COORDS(2)*pny + pny + 1
        PSZ = NODE_COORDS(3)*pnz
        PEZ = NODE_COORDS(3)*pnz + pnz + 1

        !Pick up lost points (plus ghost) on the end
        if (NODE_COORDS(1) .eq. NODE_DIMS(1)-1) then
            PEX = NX + 1
        end if
        if (NODE_COORDS(2) .eq. NODE_DIMS(2)-1) then
            PEY = NY + 1
        end if
        if (NODE_COORDS(3) .eq. NODE_DIMS(3)-1) then
            PEZ = NZ + 1
        end if
        
        !write(6,'(a,i4,a,i3,a,i3,a,i3,a)') "Rank: ", rank, " has coords: (", NODE_COORDS(1),",", &
        !    NODE_COORDS(2),",", NODE_COORDS(3),")"
        
        write(6,'(a,i4,a,i3,a,i3,a)') "Rank: ", rank, " has PSX: ", PSX,", PEX: ",PEX,"."
        write(6,'(a,i4,a,i3,a,i3,a)') "Rank: ", rank, " has PSY: ", PSY,", PEY: ",PEY,"."
        write(6,'(a,i4,a,i3,a,i3,a)') "Rank: ", rank, " has PSZ: ", PSZ,", PEZ: ",PEZ,"."


        call MPI_BARRIER(COMM_GRID, IERR)
    end subroutine
end module