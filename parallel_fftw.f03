module parallel_FFTW
    use omp_lib
    use, intrinsic :: iso_c_binding 
    include 'fftw3-mpi.f03'
    integer :: IERR_MPI, MPI_COMM_FFTW,COMM_FFTW_RANK
    type(C_PTR) :: fftw_forward_plan,fftw_backward_plan, fftw_forward_plan_dens,fftw_backward_plan_dens
    type(C_PTR) :: C_GRID_FFTW, C_GRID_T1_FFTW
    integer(C_INT) :: IERR_FFTW
    integer(C_INTPTR_T) :: alloc_local,local_NZ,local_k_offset,local_NY,local_j_offset,C_NX,C_NY,C_NZ
    contains
    subroutine init_parallel_FFTW
        implicit none
        integer :: provided
        include 'mpif.h'
        call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,IERR_MPI)
        if(provided < MPI_THREAD_FUNNELED) then
            write(6,'(a)') "ERROR: Provided MPI threading less than required!"
            call finalize_parallel_FFTW(1)
        end if
    end subroutine

    subroutine setup_parallel_topology_FFTW
        implicit none
        include 'mpif.h'
        MPI_COMM_FFTW = MPI_COMM_WORLD
        IERR_FFTW = fftw_init_threads()
        call fftw_mpi_init
        call fftw_plan_with_nthreads(omp_get_num_threads())
    end subroutine

    subroutine setup_local_allocation_fftw(NX,NY,NZ,GRID,GRID_T1)
        implicit none
        integer,intent(in)  :: NX,NY,NZ
        complex(C_DOUBLE_COMPLEX),intent(out), pointer :: GRID(:,:,:),GRID_T1(:,:,:)
        include 'mpif.h'
        C_NX = NX
        C_NY = NY
        C_NZ = NZ
        call MPI_COMM_RANK(MPI_COMM_FFTW, COMM_FFTW_RANK, IERR_MPI)
        call MPI_BARRIER(MPI_COMM_FFTW, IERR_MPI)
        alloc_local = fftw_mpi_local_size_3d(C_NZ,C_NY,C_NX,MPI_COMM_FFTW,local_NZ,local_k_offset);
		alloc_local = fftw_mpi_local_size_3d_transposed(C_NZ,C_NY,C_NX,MPI_COMM_FFTW, &
						local_NZ,local_k_offset, local_NY, local_j_offset);

        C_GRID_FFTW = fftw_alloc_complex(alloc_local)
        call c_f_pointer(C_GRID_FFTW, GRID, [C_NX, C_NY, local_NZ])

        C_GRID_T1_FFTW = fftw_alloc_complex(alloc_local)
        call c_f_pointer(C_GRID_T1_FFTW, GRID_T1, [C_NX, C_NY, local_NZ])

        if(COMM_FFTW_RANK==0) then
            write(6, '(a,i4,a,i6)') 'FFTW memory allocation complete. Local NZ is: ', local_NZ
            write(6, '(a)') 'Measuring FFT performance and selecting fastest method...'
        end if
        fftw_forward_plan = fftw_mpi_plan_dft_3d(C_NZ,C_NY,C_NX, GRID, GRID, &
                                                 MPI_COMM_FFTW, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT)
        fftw_backward_plan = fftw_mpi_plan_dft_3d(C_NZ,C_NY,C_NX, GRID, GRID, &
                                                  MPI_COMM_FFTW, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN)
        fftw_forward_plan_dens = fftw_mpi_plan_dft_3d(C_NZ,C_NY,C_NX, GRID_T1, GRID_T1, &
                                                      MPI_COMM_FFTW, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT)
        fftw_backward_plan_dens = fftw_mpi_plan_dft_3d(C_NZ,C_NY,C_NX, GRID_T1, GRID_T1, &
                                                       MPI_COMM_FFTW, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN)
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
        call fftw_free(C_GRID_T1_FFTW)
        call MPI_FINALIZE(IERR_MPI)
        if(present(status))then
            call exit(status)
        else
            call exit(0)
        end if
    end subroutine
end module
