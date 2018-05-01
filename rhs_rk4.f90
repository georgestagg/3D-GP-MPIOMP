module rhs_RK4
	use workspace
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
		
		TMPWS(2) = WS + 0.5d0*DT*TMPWS(1) !RK4_TMP2(i,j,k) = gt(i,j,k) + 0.5d0*DT*RK4_TMP1(i,j,k)
		TMPWS(3) = WS + DT*TMPWS(1)/6.0d0 !gt_out(i,j,k) = gt(i,j,k) + DT*RK4_TMP1(i,j,k)/6.0d0
		
		call halo_swap_WS(TMPWS(2)) !call halo_swap(RK4_TMP2)
		call RK4_calc_WS_RHS(TMPWS(2),TMPWS(1))

		TMPWS(2) = WS + 0.5d0*DT*TMPWS(1) !RK4_TMP2(i,j,k) = gt(i,j,k) + 0.5d0*DT*RK4_TMP1%GRID(i,j,k)
		TMPWS(3) = TMPWS(3) + DT*TMPWS(1)/3.0d0 !gt_out(i,j,k) = gt_out(i,j,k) + DT*RK4_TMP1%GRID(i,j,k)/3.0d0

		call halo_swap_WS(TMPWS(2)) !call halo_swap(RK4_TMP2%GRID)
		call RK4_calc_WS_RHS(TMPWS(2),TMPWS(1))!call RK4_gperhs(RK4_TMP2,RK4_TMP1,rt)

		TMPWS(2) = WS + DT*TMPWS(1) !RK4_TMP2(i,j,k) = gt(i,j,k) + DT*RK4_TMP1(i,j,k)
		TMPWS(3) = TMPWS(3) + DT*TMPWS(1)/3.0d0 !gt_out(i,j,k) = gt_out(i,j,k) + DT*RK4_TMP1(i,j,k)/3.0d0

		call halo_swap_WS(TMPWS(2)) !call halo_swap(RK4_TMP2)
		call RK4_calc_WS_RHS(TMPWS(2),TMPWS(1))!call RK4_gperhs(RK4_TMP2,RK4_TMP1,rt)

		TMPWS(3) = TMPWS(3) + DT*TMPWS(1)/6.0d0 !gt_out(i,j,k) = gt_out(i,j,k) + DT*RK4_TMP1(i,j,k)/6.0d0
		!All done! Copy back to main memory
		WS = TMPWS(3)

		!Renormalise WF after imaginary time decay
		if(RT .eq. 0) then
			do f = 1,FLUIDS
				call RK4_renormalise(WS%FLUID(f)%GRID,TMPWS(1)%FLUID(f)%GRID)
			end do
		end if
	end subroutine

	subroutine RK4_calc_WS_RHS(ws_in,ws_out)
		implicit none
		type(workspace_t), intent(in) :: ws_in
		type(workspace_t), intent(out) :: ws_out
		integer :: f,m
		do f = 1,FLUIDS
			call RK4_calc_RHS_Fluid(f,ws_in%FLUID(f),ws_out%FLUID(f))
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
					field_out%GRID(i,j,k) = 0.0d0
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine
	subroutine RK4_calc_RHS_Fluid(f,field_in,field_out)
		implicit none
		integer :: f,i,j,k
		class(fluid_field), intent(in) :: field_in, field_out

		if (RHSType .eq. 0) then
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						field_out%GRID(i,j,k) = -0.5d0*laplacian(field_in,i,j,k)&
							+ field_in%GRID(i,j,k)*field_in%GRID(i,j,k)*CONJG(field_in%GRID(i,j,k))&
							+ POT(i,j,k)*field_in%GRID(i,j,k) - field_in%GRID(i,j,k)&
							+ VELX*EYE*ddx(field_in,i,j,k)&
							+ VELY*EYE*ddy(field_in,i,j,k)&
							+ VELZ*EYE*ddz(field_in,i,j,k)&
							+ OMEGA*EYE*(GX(i)*ddy(field_in,i,j,k)-GY(j)*ddx(field_in,i,j,k))
					end do
				end do
			end do
			!$OMP end parallel do
		end if
		if (RHSType .eq. 1) then
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						field_out%GRID(i,j,k) = -0.5d0*laplacian(field_in,i,j,k)&
							+ harm_osc_C*field_in%GRID(i,j,k)*field_in%GRID(i,j,k)*CONJG(field_in%GRID(i,j,k))&
							+ POT(i,j,k)*field_in%GRID(i,j,k) - harm_osc_mu*field_in%GRID(i,j,k)&
							+ OMEGA*EYE*(GX(i)*ddy(field_in,i,j,k)-GY(j)*ddx(field_in,i,j,k))
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
						TMPWS(4)%FLUID(1)%GRID(i,j,k) = CONJG(WS%FLUID(f)%QP_E(i,j,k,1))*field_in%GRID(i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						field_out%GRID(i,j,k)=field_in%GRID(i,j,k)*field_in%GRID(i,j,k)*CONJG(field_in%GRID(i,j,k))&
							-field_in%GRID(i,j,k)-WS%FLUID(f)%QP_E(i,j,k,1)*d2dx2(TMPWS(4)%FLUID(1),i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						TMPWS(4)%FLUID(1)%GRID(i,j,k) = CONJG(WS%FLUID(f)%QP_E(i,j,k,2))*field_in%GRID(i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						field_out%GRID(i,j,k) = field_out%GRID(i,j,k)-WS%FLUID(f)%QP_E(i,j,k,2)&
							*d2dy2(TMPWS(4)%FLUID(1),i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz,ez
				do j = sy,ey
					do i = sx,ex
						TMPWS(4)%FLUID(1)%GRID(i,j,k) = CONJG(WS%FLUID(f)%QP_E(i,j,k,3))*field_in%GRID(i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
			!$OMP parallel do private (i,j,k) collapse(3)
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						field_out%GRID(i,j,k) = field_out%GRID(i,j,k)-WS%FLUID(f)%QP_E(i,j,k,3)&
							*d2dz2(TMPWS(4)%FLUID(1),i,j,k)
					end do
				end do
			end do
			!$OMP end parallel do
		end if
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz+1,ez-1
			do j = sy+1,ey-1
				do i = sx+1,ex-1
					field_out%GRID(i,j,k)=field_out%GRID(i,j,k)/(EYE-GAMMAC)
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine halo_swap_WS(ws_in)
		implicit none
		type(workspace_t) :: ws_in
		integer :: f,m
		include 'mpif.h'
		do f = 1,FLUIDS
			call halo_swap(ws_in%FLUID(f)%GRID,MPI_DOUBLE_COMPLEX)
		end do
		if (INC_MAG_FIELDS) then
			do m = 1,3
				call halo_swap(ws_in%MAGNETIC(f)%GRID,MPI_DOUBLE)
			end do
		end if
	end subroutine

	subroutine halo_swap(kk,mpi_type)
		use parallel
		implicit none
		integer :: mpi_type
		class(*), dimension(sx:ex,sy:ey,sz:ez) :: kk
		include 'mpif.h'
		!X dim - tag 0 and 1
		call MPI_sendrecv(kk(ex-1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), mpi_type, cart_shift(0)%rd,0,&
			kk(sx,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),mpi_type, cart_shift(0)%rs, 0, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx+1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), mpi_type, cart_shift(0)%rs,1,&
			kk(ex,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),mpi_type, cart_shift(0)%rd, 1, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		!Y dim - tag 2 and 3
		call MPI_sendrecv(kk(sx:ex,ey-1,sz:ez),(ex-sx+1)*(ez-sz+1), mpi_type, cart_shift(1)%rd,2,&
			kk(sx:ex,sy,sz:ez),(ex-sx+1)*(ez-sz+1),mpi_type, cart_shift(1)%rs, 2, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx:ex,sy+1,sz:ez),(ex-sx+1)*(ez-sz+1), mpi_type, cart_shift(1)%rs,3,&
			kk(sx:ex,ey,sz:ez),(ex-sx+1)*(ez-sz+1),mpi_type, cart_shift(1)%rd, 3, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		!Z dim - tag 4 and 5
		call MPI_sendrecv(kk(sx:ex,sy:ey,ez-1),(ex-sx+1)*(ey-sy+1), mpi_type, cart_shift(2)%rd,4,&
			kk(sx:ex,sy:ey,sz),(ex-sx+1)*(ey-sy+1),mpi_type, cart_shift(2)%rs, 4, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
		call MPI_sendrecv(kk(sx:ex,sy:ey,sz+1),(ex-sx+1)*(ey-sy+1), mpi_type, cart_shift(2)%rs,5,&
			kk(sx:ex,sy:ey,ez),(ex-sx+1)*(ey-sy+1),mpi_type, cart_shift(2)%rd, 5, MPI_COMM_GRID,MPISTAT,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
	end subroutine


	COMPLEX*16 function ddx(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(computational_field) :: field_in
		ddx = (BC(field_in,i+1,j,k)-BC(field_in,i-1,j,k))/(2.0d0*DSPACE)
	end function

	COMPLEX*16 function ddy(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(computational_field) :: field_in
		ddy = (BC(field_in,i,j+1,k)-BC(field_in,i,j-1,k))/(2.0d0*DSPACE)
	end function

	COMPLEX*16 function ddz(field_in,i,j,k)
		use params
		implicit none
		class(computational_field) :: field_in
		integer :: i,j,k
		ddz = (BC(field_in,i,j,k+1)-BC(field_in,i,j,k-1))/(2.0d0*DSPACE)
	end function

	COMPLEX*16 function d2dx2(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(computational_field) :: field_in
		select type(field_in)
		class is (magnetic_field)
			d2dx2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k))/(DSPACE**2.0d0)
		class is (fluid_field)
			d2dx2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k))/(DSPACE**2.0d0)
		end select
	end function

	COMPLEX*16 function d2dy2(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(computational_field) :: field_in
		select type(field_in)
		class is (magnetic_field)
			d2dy2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k))/(DSPACE**2.0d0)
		class is (fluid_field)
			d2dy2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k))/(DSPACE**2.0d0)
		end select
	end function

	COMPLEX*16 function d2dz2(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(computational_field) :: field_in
		select type(field_in)
		class is (magnetic_field)
			d2dz2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
		class is (fluid_field)
			d2dz2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
		end select
	end function

	COMPLEX*16 function laplacian(field_in,i,j,k)
		use params
		implicit none
		class(computational_field) :: field_in
		integer :: i,j,k
		select type(field_in)
		class is (magnetic_field)
			laplacian = (-6.0d0*field_in%GRID(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k)&
									  + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k)&
									  + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
		class is (fluid_field)
			laplacian = (-6.0d0*field_in%GRID(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k)&
									  + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k)&
									  + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
		end select
	end function

	complex*16 function BC(field_in,i,j,k)
		use params
		implicit none
		class(computational_field) :: field_in
		integer :: i,j,k,ii,jj,kk
		!Note - Ghost points means that periodic is the default

		!Zero
		if((i>=NX .or. i<=1) .and. BCX==2) then
		  BC = 0.0d0
		  RETURN
		else if ((j>=NY .or. j<=1) .and. BCY==2) then
		  BC = 0.0d0
		  RETURN
		else if ((k>=NZ .or. k<=1) .and. BCZ==2) then
		  BC = 0.0d0
		  RETURN
		end if

		!Reflective
		ii = i
		jj = j
		kk = k
		if(i == NX+1 .and. BCX == 0) ii=NX
		if(i == 0   .and. BCX == 0) ii=1
		if(j == NY+1 .and. BCY == 0) jj=NY
		if(j == 0   .and. BCY == 0) jj=1
		if(k == NZ+1 .and. BCZ == 0) kk=NZ
		if(k == 0   .and. BCZ == 0) kk=1

		select type(field_in)
		class is (magnetic_field)
			BC = field_in%GRID(ii,jj,kk)
		class is (fluid_field)
			BC = field_in%GRID(ii,jj,kk)
			!Quasi-periodic fluid field
			if(i == NX+1 .and. BCX==3) then
			  BC = BC*exp(EYE*PI*NVORTZ*GY(jj)/((NY-1)*DSPACE)-EYE*PI*NVORTY*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(i == 0 .and. BCX==3) then
			  BC =BC/exp(EYE*PI*NVORTZ*GY(jj)/((NY-1)*DSPACE)-EYE*PI*NVORTY*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(j == NY+1 .and. BCY==3) then
			  BC = BC*exp(-EYE*PI*NVORTZ*GX(ii)/((NX-1)*DSPACE)+EYE*PI*NVORTX*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(j == 0 .and. BCY==3) then
			  BC = BC/exp(-EYE*PI*NVORTZ*GX(ii)/((NX-1)*DSPACE)+EYE*PI*NVORTX*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(k == NZ+1 .and. BCZ==3) then
			  BC = BC*exp(EYE*PI*NVORTY*GX(ii)/((NX-1)*DSPACE)-EYE*PI*NVORTX*GY(jj)/((NY-1)*DSPACE)&
										+EYE*PI*NVORTX*NVORTY)
			end if
			if(k == 0 .and. BCZ==3) then
			  BC = BC/exp(EYE*PI*NVORTY*GX(ii)/((NX-1)*DSPACE)-EYE*PI*NVORTX*GY(jj)/((NY-1)*DSPACE)&
										+EYE*PI*NVORTX*NVORTY)
			end if
		end select
	end function

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