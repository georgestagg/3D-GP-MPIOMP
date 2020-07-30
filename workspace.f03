module workspace
	use params
	use parallel
	!Dynamic number of computational grids
	type computational_field
		complex(C_DOUBLE_COMPLEX), pointer :: GRID(:,:,:)
		integer :: field_number
	end type
	type, extends(computational_field) :: fluid_field
		complex(C_DOUBLE_COMPLEX), pointer :: QP_E(:,:,:,:)!Quasi-periodic variables
		complex(C_DOUBLE_COMPLEX), pointer :: DDI_K(:,:,:)!Dipole-dipole interaction variables
	end type
	type, extends(computational_field) :: magnetic_field
	end type

	type workspace_t
		type(fluid_field), pointer :: FLUID(:)
		type(magnetic_field), pointer :: MAGNETIC(:)
	end type

	type cart_shift_t
		integer :: rs,rd
	end type

	type(workspace_t) :: WS
	type(workspace_t), pointer :: TMPWS(:)

	!System-wide globals
	type(cart_shift_t),pointer :: cart_shift(:)
	real(C_DOUBLE), dimension(:,:,:), allocatable :: POT
	real(C_DOUBLE), dimension(:), allocatable :: GX,GY,GZ,KX,KY,KZ
	integer :: sx,sy,sz,ex,ey,ez,cur_step,RT
	integer :: ksx,ksy,ksz,kex,key,kez
	double precision :: TIME,DKSPACE,EDD_T1,EDD_T2
	contains
	subroutine init_workspace
		implicit none
		integer :: f,m
		if(RHSType == 0 .or. RHSType == 1) then
			call calc_local_idx_3DWithGhost(NX,NY,NZ,NGHOST,sx,ex,sy,ey,sz,ez)
			call setupCartGrid
			if(RANK .eq. 0) then
				write(6, '(a)') "Allocating distributed workspace data"
			end if
			ALLOCATE(TMPWS(3))

			ALLOCATE(WS%FLUID(FLUIDS))
			ALLOCATE(TMPWS(1)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(2)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(3)%FLUID(FLUIDS))

			do f = 1,FLUIDS
				call allocateCompGrid(WS%FLUID(f))
				call allocateCompGrid(TMPWS(1)%FLUID(f))
				call allocateCompGrid(TMPWS(2)%FLUID(f))
				call allocateCompGrid(TMPWS(3)%FLUID(f))

				WS%FLUID(f)%field_number = f
				TMPWS(1)%FLUID(f)%field_number = f
				TMPWS(2)%FLUID(f)%field_number = f
				TMPWS(3)%FLUID(f)%field_number = f
			end do
		else if(RHSType == 2) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Allocating distributed FFTW workspace data"
			end if
			ALLOCATE(TMPWS(1))

			ALLOCATE(WS%FLUID(1))
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
			call calc_local_idx_3DWithGhost(NX,NY,NZ,NGHOST,sx,ex,sy,ey,sz,ez)
			call setupCartGrid
			if(RANK .eq. 0) then
				write(6, '(a)') "Allocating distributed fluid field data"
			end if
			ALLOCATE(TMPWS(4))

			ALLOCATE(WS%FLUID(FLUIDS))
			ALLOCATE(TMPWS(1)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(2)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(3)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(4)%FLUID(1))

			do f = 1,FLUIDS
				call allocateCompGrid(WS%FLUID(f))
				call allocateCompGrid(TMPWS(1)%FLUID(f))
				call allocateCompGrid(TMPWS(2)%FLUID(f))
				call allocateCompGrid(TMPWS(3)%FLUID(f))
				WS%FLUID(f)%field_number = f
				TMPWS(1)%FLUID(f)%field_number = f
				TMPWS(2)%FLUID(f)%field_number = f
				TMPWS(3)%FLUID(f)%field_number = f
			end do
			call allocateCompGrid(TMPWS(4)%FLUID(1)) !tmp variable for quasi-periodic calcs - shared over fluids
		else if(RHSType == 4) then
			call calc_local_idx_3DWithGhost(NX,NY,NZ,NGHOST,sx,ex,sy,ey,sz,ez)
			call setupCartGrid
			if(RANK .eq. 0) then
				write(6, '(a)') "Allocating distributed fluid/magnetic workspace data"
			end if
			ALLOCATE(TMPWS(4))

			ALLOCATE(WS%FLUID(FLUIDS))
			ALLOCATE(TMPWS(1)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(2)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(3)%FLUID(FLUIDS))
			ALLOCATE(TMPWS(4)%FLUID(3))

			ALLOCATE(WS%MAGNETIC(3))
			ALLOCATE(TMPWS(1)%MAGNETIC(3))
			ALLOCATE(TMPWS(2)%MAGNETIC(3))
			ALLOCATE(TMPWS(3)%MAGNETIC(3))

			do f = 1,FLUIDS
				call allocateCompGrid(WS%FLUID(f))
				call allocateCompGrid(TMPWS(1)%FLUID(f))
				call allocateCompGrid(TMPWS(2)%FLUID(f))
				call allocateCompGrid(TMPWS(3)%FLUID(f))
				WS%FLUID(f)%field_number = f
				TMPWS(1)%FLUID(f)%field_number = f
				TMPWS(2)%FLUID(f)%field_number = f
				TMPWS(3)%FLUID(f)%field_number = f
			end do
			call allocateCompGrid(TMPWS(4)%FLUID(1)) !tmp variable for quasi-periodic calcs - shared over fluids
			INC_MAG_FIELDS = .true.
			do m = 1,3
				call allocateCompGrid(WS%MAGNETIC(m))
				call allocateCompGrid(TMPWS(1)%MAGNETIC(m))
				call allocateCompGrid(TMPWS(2)%MAGNETIC(m))
				call allocateCompGrid(TMPWS(3)%MAGNETIC(m))
				WS%MAGNETIC(m)%field_number = FLUIDS + m
				TMPWS(1)%MAGNETIC(m)%field_number = FLUIDS + m
				TMPWS(2)%MAGNETIC(m)%field_number = FLUIDS + m
				TMPWS(3)%MAGNETIC(m)%field_number = FLUIDS + m
			end do
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

		if(RHSType == 3 .or. RHSType == 4) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Pre-calculating parts of the quasi-periodic field equations"
			end if
			do f = 1,FLUIDS
				if(RANK .eq. 0) then
					write(6, '(a,i3,a,f6.3,a,f6.3,a,f6.3)') "FIELD ",f,&
						" - NVORTX is: ", NVORTALL(1,f),", NVORTY is: ", NVORTALL(2,f),", NVORTZ is: ", NVORTALL(3,f)
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
		ALLOCATE(field%GRID(sx:ex,sy:ey,sz:ez))
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
		if(RANK .eq. 0) then
			write(6, '(a)') "Pre-calculating MPI neighbours"
		end if
		ALLOCATE(cart_shift(0:2))
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
					field%QP_E(i,j,k,1) = PI*NVORTALL(2,field%field_number)*GX(i)*GZ(k)/((NX-1)*DSPACE*(NZ-1)*DSPACE)&
								  -PI*NVORTALL(3,field%field_number)*GX(i)*GY(j)/((NX-1)*DSPACE*(NY-1)*DSPACE)
				end do
			end do
		end do
		!$OMP end parallel do
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					field%QP_E(i,j,k,2) = -PI*NVORTALL(1,field%field_number)*GY(j)*GZ(k)/((NY-1)*DSPACE*(NZ-1)*DSPACE)&
								  +PI*NVORTALL(3,field%field_number)*GY(j)*GX(i)/((NY-1)*DSPACE*(NX-1)*DSPACE)
				end do
			end do
		end do
		!$OMP end parallel do
		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					field%QP_E(i,j,k,3) = PI*NVORTALL(1,field%field_number)*GZ(k)*GY(j)/((NZ-1)*DSPACE*(NY-1)*DSPACE)&
								  -PI*NVORTALL(2,field%field_number)*GZ(k)*GX(i)/((NZ-1)*DSPACE*(NX-1)*DSPACE)
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine setupGXYZ
		implicit none
		integer :: i, j ,k
		!Note: includes implicit MPI periodicity
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
end module
