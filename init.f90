module init
	use workspace
	use io
	contains
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
		else if (initialCondType .eq. 4) then
			do f = 1,FLUIDS
				call makeRandomDens(WS%FLUID(f),0.05d0)
			end do
		else
			do f = 1,FLUIDS
				WS%FLUID(f)%GRID = 1.0d0
			end do
		end if
	end subroutine

	subroutine makeRandomPhase(field)
		implicit none
		type(fluid_field) :: field
		double precision :: r1, r2
		integer :: i, j, k, n
		integer, dimension(:), allocatable :: seed

		call RANDOM_SEED(size = n)
		allocate(seed(n))
		seed = RANK
		call RANDOM_SEED(PUT = seed)

		if(RANK .eq. 0) then
			write (6, *) ' Imposing a random phase IC...'
		end if

		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					call random_number(r1)
					call random_number(r2)
					field%GRID(i,j,k) = r1*exp(2.0*PI*EYE*r2)
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine makeRandomDens(field,s)
		implicit none
		type(fluid_field) :: field
		double precision :: r1, s
		integer :: i, j, k, n
		integer, dimension(:), allocatable :: seed

		call RANDOM_SEED(size = n)
		allocate(seed(n))
		seed = RANK
		call RANDOM_SEED(PUT = seed)

		if(RANK .eq. 0) then
			write (6, *) ' Imposing a random density perturbation...'
		end if

		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					call random_number(r1)
					field%GRID(i,j,k) = 1.0d0 + r1*s - s/2.0d0
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine makeNonEquibPhase
		implicit none
		double precision :: rpKC, rpAMP, phi, ev
		integer :: i, j, k, n
		integer, dimension(:), allocatable :: seed

		call RANDOM_SEED(size = n)
		allocate(seed(n))
		seed = RANK
		call RANDOM_SEED(PUT = seed)

		ev = ((0.5*NX/PI)**2.0)*ENERV
		rpAMP = sqrt(((0.11095279993106999d0*NV**2.5d0)*(NX*NY*NZ)) / (ev**1.5d0))
		rpKC = (1.666666666666666d0*ev)/NV

		if(RANK .eq. 0) then
			write (6, *) ' Imposing the highly non-equilibrium state...'
			write (6, *) ' K-space amplitude = ', rpAMP
			write (6, *) ' Maximum wavenumber = ', sqrt(rpKC)*DKSPACE
		end if

	   !$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					if (KX(i)**2+KY(j)**2+KZ(k)**2 <= rpKC*DKSPACE**2) then
						call random_number(phi)
						WS%FLUID(1)%GRID(i,j,k) = rpAMP*exp(2.0*PI*EYE*phi)
					else
						WS%FLUID(1)%GRID(i,j,k) = 0.0d0
					end if
				end do
			end do
		end do
		!$OMP end parallel do

		call fftw_mpi_execute_dft(fftw_backward_plan, WS%FLUID(1)%GRID, WS%FLUID(1)%GRID)
		WS%FLUID(1)%GRID=WS%FLUID(1)%GRID/sqrt(dble(NX*NY*NZ))

	end subroutine
end module