module rhs_FFTW
    use workspace
	contains
    subroutine FFTW_step(rt)
        implicit none
        integer, intent(in) :: rt
        integer:: i,j,k
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    GRID_T1(i,j,k) = GRID(i,j,k)
                end do
            end do
        end do
        !$OMP end parallel do
        call fftw_mpi_execute_dft(fftw_forward_plan, GRID, GRID)
        call fftw_mpi_execute_dft(fftw_backward_plan, GRID, GRID)
        GRID = GRID/(NX*NY*NZ)
    end subroutine
end module