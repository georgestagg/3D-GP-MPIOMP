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
						GRID_T1(i,j,k) = GRID(i,j,k)*conjg(GRID(i,j,k))
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_forward_plan, GRID_T1, GRID_T1)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID_T1(i,j,k) = GRID_T1(i,j,k)*DDI_K(i,j,k)
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_backward_plan, GRID_T1, GRID_T1)
		GRID_T1 = GRID_T1/(NX*NY*NZ)


		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
				GRID(i,j,k) = GRID(i,j,k)*exp(-0.5d0*DT*EYE*(EDD_T1*GRID(i,j,k)*conjg(GRID(i,j,k)) + EDD_T2*GRID_T1(i,j,k)))
				end do
			end do
		end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_forward_plan, GRID, GRID)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID(i,j,k) = GRID(i,j,k)*exp(-0.5d0*DT*EYE*(KX(i)**2.0d0+KY(j)**2.0d0+KZ(k)**2.0d0))
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_backward_plan, GRID, GRID)
		GRID = GRID/(NX*NY*NZ)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID_T1(i,j,k) = GRID(i,j,k)*conjg(GRID(i,j,k))
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_forward_plan, GRID_T1, GRID_T1)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID_T1(i,j,k) = GRID_T1(i,j,k)*DDI_K(i,j,k)
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_backward_plan, GRID_T1, GRID_T1)
		GRID_T1 = GRID_T1/(NX*NY*NZ)
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
				GRID(i,j,k) = GRID(i,j,k)*exp(-0.5d0*DT*EYE*(EDD_T1*GRID(i,j,k)*conjg(GRID(i,j,k)) + EDD_T2*GRID_T1(i,j,k)))
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine
end module