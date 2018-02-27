module rhs_FFTW
	use workspace
	contains
	subroutine run_checks_FFTW
		implicit none
		if(BCX == 0 .or. BCY == 0 .or. BCZ == 0) then
			if(RANK .eq. 0) then
            	write(6,'(a)') "Warning: Reflective boundary conditions are not supported using the split-step fourier method."
            	write(6,'(a)') "Forcing periodic boundary conditions..."
            	BCX = 1
            	BCY = 1
            	BCZ = 1
            end if
		end if
		if(.not. (NX == NY .and. NY==NZ)) then
			if(RANK .eq. 0) then
            	write(6,'(a)') "Error: Only NX=NY=NZ is currently implemented with the split-step fourier method."
            	call finalize_parallel
            	call exit(1)
            end if
		end if
		if(GAMMAC > 0.0d0) then
			if(RANK .eq. 0) then
				write(6,'(a)') "Warning: The dissipative GPE is not supported using the split-step fourier method."
				write(6,'(a)') "Forcing GAMMAC=0.0d0"
				GAMMAC = 0.0d0
            end if
		end if
		if(OMEGA > 0.0d0) then
			if(RANK .eq. 0) then
				write(6,'(a)') "Warning: A rotating frame is not yet supported using the split-step fourier method."
				write(6,'(a)') "Forcing OMEGA=0.0d0"
				OMEGA = 0.0d0
            end if
		end if
	end subroutine

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
		call fftw_mpi_execute_dft(fftw_forward_plan_dens, GRID_T1, GRID_T1)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID_T1(i,j,k) = GRID_T1(i,j,k)/(NX*NY*NZ)*DDI_K(i,j,k)
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_backward_plan_dens, GRID_T1, GRID_T1)

		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
				GRID(i,j,k) = GRID(i,j,k)*exp(-0.5d0*DT*EYE*(POT(i,j,k) + EDD_T1*GRID(i,j,k)*conjg(GRID(i,j,k)) + EDD_T2*GRID_T1(i,j,k)))
				end do
			end do
		end do
		!$OMP end parallel do

		call fftw_mpi_execute_dft(fftw_forward_plan, GRID, GRID)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID(i,j,k) = GRID(i,j,k)/(NX*NY*NZ)*exp(-0.5d0*DT*EYE*(KX(i)**2.0d0+KZ(k)**2.0d0+KY(j)**2.0d0))
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_backward_plan, GRID, GRID)

		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID_T1(i,j,k) = GRID(i,j,k)*conjg(GRID(i,j,k))
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_forward_plan_dens, GRID_T1, GRID_T1)
		!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						GRID_T1(i,j,k) = GRID_T1(i,j,k)/(NX*NY*NZ)*DDI_K(i,j,k)
					end do
				end do
			end do
		!$OMP end parallel do
		call fftw_mpi_execute_dft(fftw_backward_plan_dens, GRID_T1, GRID_T1)

		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
				GRID(i,j,k) = GRID(i,j,k)*exp(-0.5d0*DT*EYE*(POT(i,j,k) + EDD_T1*GRID(i,j,k)*conjg(GRID(i,j,k)) + EDD_T2*GRID_T1(i,j,k)))
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine
end module