module rhs_RK4
	use workspace
	use derivs
	contains
	subroutine run_checks_RK4
		implicit none
		if(initialCondType == 2) then
			if(RANK .eq. 0) then
				write(6,'(a)') "Error: Highly non-equilibrium initialCondType is not supported with 3D block distribution."
				call finalize_parallel
				call exit(1)
			end if
		end if
	end subroutine

	subroutine RK4_step
		implicit none
		integer :: f
		call halo_swap_WS(WS)
		call RK4_calc_WS_RHS(WS,TMPWS(1))
		
		call RK4_sub_step(0.5d0,6.0d0,WS)
		
		call halo_swap_WS(TMPWS(2))
		call RK4_calc_WS_RHS(TMPWS(2),TMPWS(1))

		call RK4_sub_step(0.5d0,3.0d0,TMPWS(3))

		call halo_swap_WS(TMPWS(2))
		call RK4_calc_WS_RHS(TMPWS(2),TMPWS(1))

		call RK4_sub_step(1.0d0,3.0d0,TMPWS(3))

		call halo_swap_WS(TMPWS(2))
		call RK4_calc_WS_RHS(TMPWS(2),TMPWS(1))

		call RK4_final_step

		!Renormalise WF after imaginary time decay
		if(RT .eq. 0) then
			do f = 1,FLUIDS
				call RK4_renormalise_fluid(WS%FLUID(f)%GRID,TMPWS(1)%FLUID(f)%GRID)
			end do
		end if
	end subroutine

	subroutine RK4_sub_step(d1,d2,IWS)
		implicit none
		integer :: m,f,i,j,k
		double precision :: d1,d2
		type(workspace_t),intent(in) :: IWS
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						TMPWS(2)%FLUID(f)%GRID(i,j,k) = WS%FLUID(f)%GRID(i,j,k) + TMPWS(1)%FLUID(f)%GRID(i,j,k)*d1*DT 
						TMPWS(3)%FLUID(f)%GRID(i,j,k) = IWS%FLUID(f)%GRID(i,j,k) + TMPWS(1)%FLUID(f)%GRID(i,j,k)*(DT/d2)
					end do
				end do
			end do
		end do
		!$OMP end parallel do
		if (INC_MAG_FIELDS) then
			!$OMP parallel do private (m,i,j,k) collapse(4)
			do m = 1,3
				do k = sz+1,ez-1
					do j = sy+1,ey-1
						do i = sx+1,ex-1
							TMPWS(2)%MAGNETIC(m)%GRID_R(i,j,k) = WS%MAGNETIC(m)%GRID_R(i,j,k) + TMPWS(1)%MAGNETIC(m)%GRID_R(i,j,k)*d1*DT 
							TMPWS(3)%MAGNETIC(m)%GRID_R(i,j,k) = IWS%MAGNETIC(m)%GRID_R(i,j,k) + TMPWS(1)%MAGNETIC(m)%GRID_R(i,j,k)*(DT/d2)
						end do
					end do
				end do
			end do
			!$OMP end parallel do
		end if
	end subroutine

	subroutine RK4_final_step
		implicit none
		integer :: m,f,i,j,k
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WS%FLUID(f)%GRID(i,j,k) = TMPWS(3)%FLUID(f)%GRID(i,j,k) + TMPWS(1)%FLUID(f)%GRID(i,j,k)*(DT/6.0d0)
					end do
				end do
			end do
		end do
		!$OMP end parallel do
		if (INC_MAG_FIELDS) then
			!$OMP parallel do private (m,i,j,k) collapse(4)
			do m = 1,3
				do k = sz+1,ez-1
					do j = sy+1,ey-1
						do i = sx+1,ex-1
							WS%MAGNETIC(m)%GRID_R(i,j,k) = TMPWS(3)%MAGNETIC(m)%GRID_R(i,j,k) + TMPWS(1)%MAGNETIC(m)%GRID_R(i,j,k)*(DT/6.0d0)
						end do
					end do
				end do
			end do
			!$OMP end parallel do
		end if
	end subroutine

	subroutine RK4_calc_WS_RHS(ws_in,ws_out)
		implicit none
		type(workspace_t), intent(in) :: ws_in
		type(workspace_t), intent(out) :: ws_out
		integer :: f,m
		do f = 1,FLUIDS
			call RK4_calc_RHS_Fluid(f,ws_in,ws_out)
		end do
		if (INC_MAG_FIELDS) then
			do m = 1,3
				call RK4_calc_RHS_Magnetic(m,ws_in%MAGNETIC(m),ws_out%MAGNETIC(m))
			end do
		end if
	end subroutine
	subroutine RK4_calc_RHS_Magnetic(m,field_in,field_out)
		implicit none
		integer :: m,i,j,k
		class(magnetic_field), intent(in) :: field_in, field_out
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz+1,ez-1
			do j = sy+1,ey-1
				do i = sx+1,ex-1
					field_out%GRID_R(i,j,k) = 0.0d0
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine
	subroutine RK4_calc_RHS_Fluid(f,ws_in,ws_out)
		implicit none
		integer :: f,i,j,k,p
		type(workspace_t), intent(in) :: ws_in
		type(workspace_t), intent(out) :: ws_out

		if (RHSType .eq. 0) then
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						ws_out%FLUID(f)%GRID(i,j,k) = -0.5d0*laplacian(ws_in%FLUID(f),i,j,k)&
							+ POT(i,j,k)*ws_in%FLUID(f)%GRID(i,j,k) - ws_in%FLUID(f)%GRID(i,j,k)&
							+ VELX*EYE*ddx(ws_in%FLUID(f),i,j,k)&
							+ VELY*EYE*ddy(ws_in%FLUID(f),i,j,k)&
							+ VELZ*EYE*ddz(ws_in%FLUID(f),i,j,k)&
							+ OMEGA*EYE*(GX(i)*ddy(ws_in%FLUID(f),i,j,k)-GY(j)*ddx(ws_in%FLUID(f),i,j,k))
					end do
				end do
			end do
			!!$OMP end parallel do
			do p = 1,FLUIDS
				!$OMP parallel do private (i,j,k) collapse(3)
				do k = sz+1,ez-1
					do j = sy+1,ey-1
						do i = sx+1,ex-1
							ws_out%FLUID(f)%GRID(i,j,k) = ws_out%FLUID(f)%GRID(i,j,k) +&
								GG(f,p)*ws_in%FLUID(p)%GRID(i,j,k)*CONJG(ws_in%FLUID(p)%GRID(i,j,k))*ws_in%FLUID(f)%GRID(i,j,k)
						end do
					end do
				end do
				!$OMP end parallel do
			end do
		end if
		if (RHSType .eq. 1) then
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						ws_out%FLUID(f)%GRID(i,j,k) = -0.5d0*laplacian(ws_in%FLUID(f),i,j,k)&
							+ harm_osc_C*ws_in%FLUID(f)%GRID(i,j,k)*ws_in%FLUID(f)%GRID(i,j,k)*CONJG(ws_in%FLUID(f)%GRID(i,j,k))&
							+ POT(i,j,k)*ws_in%FLUID(f)%GRID(i,j,k) - harm_osc_mu*ws_in%FLUID(f)%GRID(i,j,k)&
							+ OMEGA*EYE*(GX(i)*ddy(ws_in%FLUID(f),i,j,k)-GY(j)*ddx(ws_in%FLUID(f),i,j,k))
					end do
				end do
			end do
			!$OMP end parallel do
		end if
		if (RHSType .eq. 3) then
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						TMPWS(4)%FLUID(1)%GRID(i,j,k) = CONJG(WS%FLUID(f)%QP_E(i,j,k,1))*ws_in%FLUID(f)%GRID(i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						ws_out%FLUID(f)%GRID(i,j,k)= -ws_in%FLUID(f)%GRID(i,j,k)-WS%FLUID(f)%QP_E(i,j,k,1)&
							*d2dx2(TMPWS(4)%FLUID(1),i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						TMPWS(4)%FLUID(1)%GRID(i,j,k) = CONJG(WS%FLUID(f)%QP_E(i,j,k,2))*ws_in%FLUID(f)%GRID(i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						ws_out%FLUID(f)%GRID(i,j,k) = ws_out%FLUID(f)%GRID(i,j,k)-WS%FLUID(f)%QP_E(i,j,k,2)&
							*d2dy2(TMPWS(4)%FLUID(1),i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						TMPWS(4)%FLUID(1)%GRID(i,j,k) = CONJG(WS%FLUID(f)%QP_E(i,j,k,3))*ws_in%FLUID(f)%GRID(i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						ws_out%FLUID(f)%GRID(i,j,k) = ws_out%FLUID(f)%GRID(i,j,k)-WS%FLUID(f)%QP_E(i,j,k,3)&
							*d2dz2(TMPWS(4)%FLUID(1),i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			
			do p = 1,FLUIDS
				!$OMP parallel do private (i,j,k) collapse(3)
				do k = sz+1,ez-1
					do j = sy+1,ey-1
						do i = sx+1,ex-1
							ws_out%FLUID(f)%GRID(i,j,k) = ws_out%FLUID(f)%GRID(i,j,k)&
								+GG(f,p)*ws_in%FLUID(p)%GRID(i,j,k)*CONJG(ws_in%FLUID(p)%GRID(i,j,k))*ws_in%FLUID(f)%GRID(i,j,k)
						end do
					end do
				end do
				!$OMP end parallel do
			end do

		end if
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz+1,ez-1
			do j = sy+1,ey-1
				do i = sx+1,ex-1
					ws_out%FLUID(f)%GRID(i,j,k)=ws_out%FLUID(f)%GRID(i,j,k)/(EYE-GAMMAC)
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine halo_swap_WS(ws_in)
		implicit none
		type(workspace_t) :: ws_in
		integer :: f,m
		do f = 1,FLUIDS
			call halo_swap_complex(ws_in%FLUID(f)%GRID)
		end do
		if (INC_MAG_FIELDS) then
			do m = 1,3
				call halo_swap_real(ws_in%MAGNETIC(f)%GRID_R)
			end do
		end if
	end subroutine

	subroutine halo_swap_complex(kk)
		use parallel
		implicit none
		complex*16, dimension(sx:ex,sy:ey,sz:ez) :: kk
		include 'mpif.h'
		!X dim - tag 0 and 1
		call MPI_sendrecv(kk(ex-1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, cart_shift(0)%rd,0,&
			kk(sx,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, cart_shift(0)%rs, 0, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx+1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, cart_shift(0)%rs,1,&
			kk(ex,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, cart_shift(0)%rd, 1, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		!Y dim - tag 2 and 3
		call MPI_sendrecv(kk(sx:ex,ey-1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, cart_shift(1)%rd,2,&
			kk(sx:ex,sy,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, cart_shift(1)%rs, 2, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx:ex,sy+1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, cart_shift(1)%rs,3,&
			kk(sx:ex,ey,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, cart_shift(1)%rd, 3, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		!Z dim - tag 4 and 5
		call MPI_sendrecv(kk(sx:ex,sy:ey,ez-1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE_COMPLEX, cart_shift(2)%rd,4,&
			kk(sx:ex,sy:ey,sz),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE_COMPLEX, cart_shift(2)%rs, 4, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx:ex,sy:ey,sz+1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE_COMPLEX, cart_shift(2)%rs,5,&
			kk(sx:ex,sy:ey,ez),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE_COMPLEX, cart_shift(2)%rd, 5, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
	end subroutine

	subroutine halo_swap_real(kk)
		use parallel
		implicit none
		double precision, dimension(sx:ex,sy:ey,sz:ez) :: kk
		include 'mpif.h'
		!X dim - tag 0 and 1
		call MPI_sendrecv(kk(ex-1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE, cart_shift(0)%rd,0,&
			kk(sx,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE, cart_shift(0)%rs, 0, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx+1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE, cart_shift(0)%rs,1,&
			kk(ex,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE, cart_shift(0)%rd, 1, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		!Y dim - tag 2 and 3
		call MPI_sendrecv(kk(sx:ex,ey-1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE, cart_shift(1)%rd,2,&
			kk(sx:ex,sy,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE, cart_shift(1)%rs, 2, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx:ex,sy+1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE, cart_shift(1)%rs,3,&
			kk(sx:ex,ey,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE, cart_shift(1)%rd, 3, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		!Z dim - tag 4 and 5
		call MPI_sendrecv(kk(sx:ex,sy:ey,ez-1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE, cart_shift(2)%rd,4,&
			kk(sx:ex,sy:ey,sz),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE, cart_shift(2)%rs, 4, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx:ex,sy:ey,sz+1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE, cart_shift(2)%rs,5,&
			kk(sx:ex,sy:ey,ez),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE, cart_shift(2)%rd, 5, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
	end subroutine

	subroutine RK4_renormalise_fluid(gt,tgt)
		implicit none
		complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt,tgt
		include 'mpif.h'
		double precision :: local_norm,total_norm,local_pot_int,total_pot_int

		tgt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1) = gt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1)*&
												conjg(gt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))

		local_norm=sum(tgt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))*DSPACE*DSPACE*DSPACE
		call MPI_Allreduce(local_norm, total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW,IERR)

		tgt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1) = 0.0d0
		where (POT<1.0d0) tgt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1)=1.0d0

		local_pot_int=sum(tgt(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))*DSPACE*DSPACE*DSPACE
		call MPI_Allreduce(local_pot_int, total_pot_int, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW,IERR)

		gt(sx:ex,sy:ey,sz:ez) = gt(sx:ex,sy:ey,sz:ez)/sqrt(total_norm)*sqrt(total_pot_int)
	end subroutine
end module