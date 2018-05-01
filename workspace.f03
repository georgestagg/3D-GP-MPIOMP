module workspace
	use params
	use parallel

	!Dynamic number of computational grids
	type computational_field
	end type
	type, extends(computational_field) :: fluid_field
		complex(C_DOUBLE_COMPLEX), pointer :: GRID(:,:,:) !Fluid wavefunction
		complex(C_DOUBLE_COMPLEX), pointer :: QP_E(:,:,:,:) !Quasi-periodic variables
		complex(C_DOUBLE_COMPLEX), pointer :: DDI_K(:,:,:) !Dipole-dipole interaction variables
	end type
	type, extends(computational_field) :: magnetic_field
		real(C_DOUBLE), pointer :: GRID(:,:,:)
	end type

	type workspace_t
		type(fluid_field), pointer :: FLUID(:)
		type(magnetic_field), pointer :: MAGNETIC(:)
	end type

	type cart_shift_t
		integer :: rs,rd
	end type

	interface operator (+)
		procedure workspace_add_WS_WS
	end interface

	interface operator (*)
		procedure workspace_mult_WS_WS
		procedure workspace_mult_WS_SCALAR
		procedure workspace_mult_WS_SCALARC
		procedure workspace_mult_SCALAR_WS
		procedure workspace_mult_SCALARC_WS
	end interface

	interface operator (/)
		procedure workspace_div_WS_SCALAR
		procedure workspace_div_WS_SCALARC
		procedure workspace_div_SCALAR_WS
		procedure workspace_div_SCALARC_WS
	end interface

	type(workspace_t) :: WS
	type(workspace_t), pointer :: TMPWS(:)

	!System-wide globals
	type(cart_shift_t),pointer :: cart_shift(:)
	real(C_DOUBLE), dimension(:,:,:), allocatable :: POT
	real(C_DOUBLE), dimension(:,:), allocatable :: GG
	real(C_DOUBLE), dimension(:), allocatable :: GX,GY,GZ,KX,KY,KZ
	integer :: sx,sy,sz,ex,ey,ez,cur_step,RT
	double precision :: TIME,DKSPACE,EDD_T1,EDD_T2
	contains
	subroutine init_workspace
		implicit none
		integer :: f,m
		if(RANK .eq. 0) then
			write(6, '(a)') "Allocating distributed workspace data"
		end if
		if(RHSType == 0 .or. RHSType == 1) then
			call calc_local_idx_3DWithGhost(NX,NY,NZ,sx,ex,sy,ey,sz,ez)
			ALLOCATE(cart_shift(0:2))
			call setupCartGrid
			ALLOCATE(TMPWS(3))
			ALLOCATE(WS%FLUID(FLUIDS))
			ALLOCATE(GG(FLUIDS,FLUIDS))
			do f = 1,FLUIDS
				call allocateCompGri(WS%FLUID(f))
				call allocateCompGri(TMPWS(1)%FLUID(f))
				call allocateCompGri(TMPWS(2)%FLUID(f))
				call allocateCompGri(TMPWS(3)%FLUID(f))
			end do
		else if(RHSType == 2) then
			ALLOCATE(WS%FLUID(1))
			ALLOCATE(TMPWS(1))
			ALLOCATE(TMPWS(1)%FLUID(1))
			call setup_local_allocation_fftw(NX,NY,NZ,WS%FLUID(1)%GRID,TMPWS(1)%FLUID(1)%GRID)
			call parallel_barrier
			sx = 1
			ex = NX
			sy = 1
			ey = NY
			sz = 1
			ez = local_NZ
		else if(RHSType == 3) then
			call calc_local_idx_3DWithGhost(NX,NY,NZ,sx,ex,sy,ey,sz,ez)
			ALLOCATE(cart_shift(0:2))
			call setupCartGrid
			ALLOCATE(WS%FLUID(FLUIDS))
			ALLOCATE(GG(FLUIDS,FLUIDS))
			ALLOCATE(TMPWS(4))
			do f = 1,FLUIDS
				call allocateCompGrid(WS%FLUID(f))
				call allocateCompGri(TMPWS(1)%FLUID(f))
				call allocateCompGri(TMPWS(2)%FLUID(f))
				call allocateCompGri(TMPWS(3)%FLUID(f))
			end do
			call allocateCompGri(TMPWS(4)%FLUID(1)) !tmp variable for quasi-periodic calcs - shared over fluids
			if (INC_MAG_FIELDS) then
				ALLOCATE(WS%MAGNETIC(3))
				do m = 1,3
					call allocateCompGrid(WS%MAGNETIC(m))
					call allocateCompGri(TMPWS(1)%MAGNETIC(m))
					call allocateCompGri(TMPWS(2)%MAGNETIC(m))
					call allocateCompGri(TMPWS(3)%MAGNETIC(m))
				end do
			end if
		end if

		if(RANK .eq. 0) then
			write(6, '(a)') "Calculating distributed position and momentum space coordinates"
		end if
		ALLOCATE(POT(sx:ex,sy:ey,sz:ez))
		ALLOCATE(GX(sx:ex))
		ALLOCATE(GY(sy:ey))
		ALLOCATE(GZ(sz:ez))
		ALLOCATE(KX(sx:ex))
		ALLOCATE(KY(sy:ey))
		ALLOCATE(KZ(sz:ez))
		call setupGXYZ
		call setupKXYZ
		call parallel_barrier

		if(RHSType == 2) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Pre-calculating dipole-dipole k-space potential"
				write(6, '(a,f6.3)') "EDD is: ", EDD
			end if
			EDD_T1 = 1.0d0/(1.0d0-EDD)
			EDD_T2 = EDD/(1.0d0-EDD)
			ALLOCATE(WS%FLUID(1)%DDI_K(sx:ex,sy:ey,sz:ez))
			call setupDDIK(WS%FLUID(1))
		end if

		if(RHSType == 3) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Pre-calculating parts of the quasi-periodic field equations"
			end if
			do f = 1,FLUIDS
				if(RANK .eq. 0) then
					write(6, '(a,i3,a,f6.3,a,f6.3,a,f6.3)') "FIELD ",f,&
						" - NVORTX is: ", NVORTX,", NVORTY is: ", NVORTY,", NVORTZ is: ", NVORTZ
				end if
				ALLOCATE(WS%FLUID(f)%QP_E(sx:ex,sy:ey,sz:ez,3))
				call setupEXEYEZ(WS%FLUID(f))
			end do
		end if
		call parallel_barrier
	end subroutine

	subroutine allocateCompGrid(field)
		implicit none
		class(computational_field) :: field
		select type(field)
		class is (fluid_field)
			ALLOCATE(field%GRID(sx:ex,sy:ey,sz:ez))
		class is (magnetic_field)
			ALLOCATE(field%GRID(sx:ex,sy:ey,sz:ez))
		end select
	end subroutine

	subroutine setupDDIK(field)
		implicit none
		integer :: i, j ,k
		type(fluid_field) :: field
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					!NOTE: This term is calculated in transposed k-space, so KY is really KZ.
					!NB: http://fftw.org/doc/Transposed-distributions.html#Transposed-distributions
					field%DDI_K(i,j,k) = 3.0d0*(KY(j)**2.0d0)/(KX(i)**2.0d0+KY(j)**2.0d0+KZ(k)**2.0d0+1e-20) - 1.0d0
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine setupCartGrid
		implicit none
		include 'mpif.h'
		call MPI_Cart_shift(MPI_COMM_GRID,0,1,cart_shift(0)%rs,cart_shift(0)%rd,IERR)
		call MPI_Cart_shift(MPI_COMM_GRID,1,1,cart_shift(1)%rs,cart_shift(1)%rd,IERR)
		call MPI_Cart_shift(MPI_COMM_GRID,2,1,cart_shift(2)%rs,cart_shift(2)%rd,IERR)
		if(IERR .ne. MPI_SUCCESS) then
			write(6,*) "Error running MPI_Cart_shift. Error code: ",IERR
			write(6,*) "Something has gone wrong... Quitting."
			CALL EXIT(1)
		end if
	end subroutine

	subroutine setupEXEYEZ(field)
		implicit none
		integer :: i, j ,k
		type(fluid_field) :: field
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					field%QP_E(i,j,k,1) = exp(EYE*PI*NVORTY*GX(i)*GZ(k)/((NX-1)*DSPACE*(NZ-1)*DSPACE)&
								  -EYE*PI*NVORTZ*GX(i)*GY(j)/((NX-1)*DSPACE*(NY-1)*DSPACE))
				end do
			end do
		end do
		!$OMP end parallel do
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					field%QP_E(i,j,k,2) = exp(-EYE*PI*NVORTX*GY(j)*GZ(k)/((NY-1)*DSPACE*(NZ-1)*DSPACE)&
								  +EYE*PI*NVORTZ*GY(j)*GX(i)/((NY-1)*DSPACE*(NX-1)*DSPACE))
				end do
			end do
		end do
		!$OMP end parallel do
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					field%QP_E(i,j,k,3) = exp(EYE*PI*NVORTX*GZ(k)*GY(j)/((NZ-1)*DSPACE*(NY-1)*DSPACE)&
								  -EYE*PI*NVORTY*GZ(k)*GX(i)/((NZ-1)*DSPACE*(NX-1)*DSPACE))
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine setupGXYZ
		implicit none
		integer :: i, j ,k
		do k = sz,ez
			GZ(k) = (k+local_k_offset-1)*DSPACE-ZSHIFT
		end do		

		do j = sy,ey
			GY(j) = (j-1)*DSPACE-YSHIFT
		end do

		do i = sx,ex
			GX(i) = (i-1)*DSPACE-XSHIFT
		end do
	end subroutine

	subroutine setupKXYZ
		implicit none
		integer :: i, j ,k
		DKSPACE = 2*PI/(DSPACE*NX)

		do k = sz,ez
			if ((k+local_k_offset-1) < NZ/2) then
				KZ(k) = 2*PI*(k+local_k_offset-1)/(DSPACE*NZ)
			else
				KZ(k) = -2*PI*(NZ-k-local_k_offset+1)/(DSPACE*NZ)
			end if
		end do

		do j = sy,ey
			if (j-1 < NY/2) then
				KY(j) = 2*PI*(j-1)/(DSPACE*NY)
			else
				KY(j) = -2*PI*(NY-j+1)/(DSPACE*NY)
			end if
		end do

		do i = sx,ex
			if (i-1 < NX/2) then
				KX(i) = 2*PI*(i-1)/(DSPACE*NX)
			else
				KX(i) = -2*PI*(NX-i+1)/(DSPACE*NX)
			end if
		end do
	end subroutine

	!The remaining functions are operator overloading of various types
	function workspace_add_WS_WS(WS1,WS2) result (WSO)
		implicit none
		integer :: f,i,j,k
		type(workspace_t), intent(in) :: WS1,WS2
		type(workspace_t) :: WSO
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WSO%FLUID(f)%GRID(i,j,k) = WS1%FLUID(f)%GRID(i,j,k)+WS2%FLUID(f)%GRID(i,j,k)
					end do
				end do
			end do
		end do
		!$OMP end parallel do
	end function

	function workspace_mult_WS_WS(WS1,WS2) result (WSO)
		implicit none
		integer :: f,i,j,k
		type(workspace_t), intent(in) :: WS1,WS2
		type(workspace_t) :: WSO
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WSO%FLUID(f)%GRID(i,j,k) = WS1%FLUID(f)%GRID(i,j,k)*WS2%FLUID(f)%GRID(i,j,k)
					end do
				end do
			end do
		end do
		!$OMP end parallel do
	end function

	function workspace_mult_WS_SCALAR(WS,SC) result (WSO)
		implicit none
		integer :: f,i,j,k
		type(workspace_t), intent(in) :: WS
		real(C_DOUBLE), intent(in) :: SC
		type(workspace_t) :: WSO
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WSO%FLUID(f)%GRID(i,j,k) = WS%FLUID(f)%GRID(i,j,k)*SC
					end do
				end do
			end do
		end do
		!$OMP end parallel do
	end function

	function workspace_mult_WS_SCALARC(WS,SCC) result (WSO)
		implicit none
		integer :: f,i,j,k
		type(workspace_t), intent(in) :: WS
		complex(C_DOUBLE_COMPLEX), intent(in) :: SCC
		type(workspace_t) :: WSO
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WSO%FLUID(f)%GRID(i,j,k) = WS%FLUID(f)%GRID(i,j,k)*SCC
					end do
				end do
			end do
		end do
		!$OMP end parallel do
	end function

	function workspace_div_WS_SCALAR(WS,SC) result (WSO)
		implicit none
		integer :: f,i,j,k
		type(workspace_t), intent(in) :: WS
		real(C_DOUBLE), intent(in) :: SC
		type(workspace_t) :: WSO
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WSO%FLUID(f)%GRID(i,j,k) = WS%FLUID(f)%GRID(i,j,k)/SC
					end do
				end do
			end do
		end do
		!$OMP end parallel do
	end function

	function workspace_div_WS_SCALARC(WS,SCC) result (WSO)
		implicit none
		integer :: f,i,j,k
		type(workspace_t), intent(in) :: WS
		complex(C_DOUBLE_COMPLEX), intent(in) :: SCC
		type(workspace_t) :: WSO
		!$OMP parallel do private (f,i,j,k) collapse(4)
		do f = 1,FLUIDS
			do k = sz+1,ez-1
				do j = sy+1,ey-1
					do i = sx+1,ex-1
						WSO%FLUID(f)%GRID(i,j,k) = WS%FLUID(f)%GRID(i,j,k)/SCC
					end do
				end do
			end do
		end do
		!$OMP end parallel do
	end function

	function workspace_mult_SCALAR_WS(SC,WS) result (WSO)
		implicit none
		type(workspace_t), intent(in) :: WS
		real(C_DOUBLE), intent(in) :: SC
		type(workspace_t) :: WSO
		WSO = workspace_mult_WS_SCALAR(WS,SC)
	end function
	function workspace_mult_SCALARC_WS(SCC,WS) result (WSO)
		implicit none
		type(workspace_t), intent(in) :: WS
		complex(C_DOUBLE_COMPLEX), intent(in) :: SCC
		type(workspace_t) :: WSO
		WSO = workspace_mult_WS_SCALARC(WS,SCC)
	end function
	function workspace_div_SCALAR_WS(SC,WS) result (WSO)
		implicit none
		type(workspace_t), intent(in) :: WS
		real(C_DOUBLE), intent(in) :: SC
		type(workspace_t) :: WSO
		WSO = workspace_div_WS_SCALAR(WS,SC)
	end function
	function workspace_div_SCALARC_WS(SCC,WS) result (WSO)
		implicit none
		type(workspace_t), intent(in) :: WS
		complex(C_DOUBLE_COMPLEX), intent(in) :: SCC
		type(workspace_t) :: WSO
		WSO = workspace_div_WS_SCALARC(WS,SCC)
	end function

subroutine initCond
	implicit none
	integer :: f
	if(RANK .eq. 0) then
		write (6, *) 'Applying initial condition...'
	end if
	TIME = 0.0d0
	if(initialCondType .eq. 0) then
		do f = 1,FLUIDS
			WS%FLUID(f)%GRID = 1.0d0
		end do
	else if (initialCondType .eq. 1) then
		do f = 1,FLUIDS
			call makeRandomPhase(WS%FLUID(f))
		end do
	else if (initialCondType .eq. 2) then
		call makeNonEquibPhase
	else if (initialCondType .eq. 3) then
		call read_wf_file(ICRfilename)
	else
		do f = 1,FLUIDS
			WS%FLUID(f)%GRID = 1.0d0
		end do
	end if
end subroutine

end module