module parallel_3DWithGhost
    use omp_lib
    integer :: NODE_COORDS(3),NODE_DIMS(3),MPI_COMM_GRID,COMM_GRID_RANK,IERR_3DWG
    integer, dimension(:), ALLOCATABLE :: MPISTAT
    contains
    subroutine init_parallel_3DWithGhost
        implicit none
        include 'mpif.h'
        ALLOCATE(MPISTAT(MPI_STATUS_SIZE))
        call MPI_INIT(IERR_3DWG)
        call MPI_BARRIER(MPI_COMM_WORLD, IERR_3DWG)
    end subroutine

    subroutine setup_parallel_topology_3DWithGhost(NPROCS,RANK)
        implicit none
        include 'mpif.h'
        integer,intent(in) :: NPROCS,RANK
        logical :: PERIODIC(3)
        !These are not for setting GPE periodicity, it is the periodicity of the MPI topology
        PERIODIC(1) = .true.
        PERIODIC(2) = .true.
        PERIODIC(3) = .true.
        NODE_DIMS = 0
        call MPI_DIMS_CREATE(NPROCS, 3, NODE_DIMS,IERR_3DWG)
        call MPI_CART_CREATE(MPI_COMM_WORLD, 3, NODE_DIMS, PERIODIC, .true., MPI_COMM_GRID, IERR_3DWG)
        if(RANK .eq. 0) then
            write(6,'(a,i6,a,i6,a,i6,a)') "Topology is: (", NODE_DIMS(1),",", &
                NODE_DIMS(2),",", NODE_DIMS(3),")."
        end if
        call MPI_COMM_RANK(MPI_COMM_GRID, COMM_GRID_RANK, IERR_3DWG)
        call MPI_BARRIER(MPI_COMM_GRID, IERR_3DWG)
    end subroutine

    !get local grid sizes
    subroutine calc_local_idx_3DWithGhost(NX,NY,NZ,PSX,PEX,PSY,PEY,PSZ,PEZ)
        implicit none
        integer :: pnx,pny,pnz
        integer,intent(in)  :: NX,NY,NZ
        integer,intent(out) :: PSX,PEX,PSY,PEY,PSZ,PEZ
        pnx = NX/NODE_DIMS(1)
        pny = NY/NODE_DIMS(2)
        pnz = NY/NODE_DIMS(3)     
        call MPI_CART_COORDS(MPI_COMM_GRID,COMM_GRID_RANK,3,NODE_COORDS,IERR_3DWG)

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
        call MPI_BARRIER(MPI_COMM_GRID, IERR_3DWG)
    end subroutine

    subroutine finalize_parallel_3DWithGhost
        implicit none
        call MPI_BARRIER(MPI_COMM_GRID, IERR_3DWG)
        call FLUSH()
        call MPI_FINALIZE(IERR_3DWG)
    end subroutine
end module