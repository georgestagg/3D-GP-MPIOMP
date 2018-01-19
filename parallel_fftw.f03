module parallel_FFTW
    use omp_lib
    use, intrinsic :: iso_c_binding 
    include 'fftw3-mpi.f03'
    integer :: IERR_MPI, MPI_COMM_FFTW,COMM_FFTW_RANK
    type(C_PTR) :: fftw_forward_plan,fftw_backward_plan, C_GRID_FFTW
    integer(C_INT) :: IERR_FFTW
    integer(C_INTPTR_T) :: alloc_local,local_NZ,local_k_offset,C_NX,C_NY,C_NZ
    contains
    subroutine init_parallel_FFTW
        implicit none
        integer :: provided
        include 'mpif.h'
        MPI_COMM_FFTW = MPI_COMM_WORLD
        call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,IERR_MPI)
        if(provided < MPI_THREAD_FUNNELED) then
            write(6,'(a)') "ERROR: Provided MPI threading less than required!"
            call finalize_parallel_FFTW(1)
        end if
    end subroutine

    subroutine setup_parallel_topology_FFTW
        implicit none

        IERR_FFTW = fftw_init_threads()
        call fftw_mpi_init
        call fftw_plan_with_nthreads(omp_get_num_threads())
    end subroutine

    subroutine setup_local_allocation(NX,NY,NZ,GRID)
        implicit none
        integer,intent(in)  :: NX,NY,NZ
        complex(C_DOUBLE_COMPLEX),intent(in), pointer :: GRID(:,:,:)
        include 'mpif.h'
        C_NX = NX
        C_NY = NY
        C_NZ = NZ
        call MPI_COMM_RANK(MPI_COMM_FFTW, COMM_FFTW_RANK, IERR_MPI)
        alloc_local = fftw_mpi_local_size_3d(C_NZ,C_NY,C_NX,MPI_COMM_FFTW,local_NZ,local_k_offset);
        call MPI_BARRIER(MPI_COMM_FFTW, IERR_MPI)
        write(6, '(a,i4,a,i6)') 'I am rank: ', COMM_FFTW_RANK,', my local NZ is: ', local_NZ
        call MPI_BARRIER(MPI_COMM_FFTW, IERR_MPI)
        if(COMM_FFTW_RANK==0) then
            write(6, '(a)') 'Measuring FFT performance and selecting fastest method...'
        end if
        call MPI_BARRIER(MPI_COMM_FFTW, IERR_MPI)
        C_GRID_FFTW = fftw_alloc_complex(alloc_local)
        call c_f_pointer(C_GRID_FFTW, GRID, [C_NX, C_NY, local_NZ])
        fftw_forward_plan = fftw_mpi_plan_dft_3d(C_NZ,C_NY,C_NX, GRID, GRID, MPI_COMM_FFTW, FFTW_FORWARD, FFTW_PATIENT)
        fftw_backward_plan = fftw_mpi_plan_dft_3d(C_NZ,C_NY,C_NX, GRID, GRID, MPI_COMM_FFTW, FFTW_BACKWARD, FFTW_PATIENT)
        if(COMM_FFTW_RANK==0) then
            write(6, '(a)') 'Done!'
        end if
        call MPI_BARRIER(MPI_COMM_FFTW, IERR_MPI)
    end subroutine

    subroutine finalize_parallel_FFTW(status)
        implicit none
        integer, optional :: status
        call fftw_destroy_plan(fftw_forward_plan)
        call fftw_destroy_plan(fftw_backward_plan)
        call fftw_mpi_cleanup
        call fftw_free(C_GRID_FFTW)
        call MPI_FINALIZE(IERR_MPI)
        if(present(status))then
            call exit(status)
        else
            call exit(0)
        end if
    end subroutine
end module