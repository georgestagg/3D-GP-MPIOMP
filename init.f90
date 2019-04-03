module init
	use workspace
	use io
	contains
	subroutine initCond
		integer :: f
		if(RANK .eq. 0) then
			write (6, "(a)") 'Applying initial condition...'
		end if
		TIME = 0.0d0
		if (initialCondType .eq. 1) then
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
		else if (initialCondType .eq. 5) then
			include 'ic.in'
		else
			do f = 1,FLUIDS
				call makeConst(WS%FLUID(f),(1.0d0,0.0d0))
			end do
		end if
	end subroutine

	subroutine makeConst(field,c)
		implicit none
		type(fluid_field) :: field
		integer :: i, j, k
		complex*16 :: c

		!$OMP parallel do private (i,j,k) collapse(3)
		do k = sz,ez
			do j = sy,ey
				do i = sx,ex
					field%GRID(i,j,k) = c
				end do
			end do
		end do
		!$OMP end parallel do
	end subroutine

	subroutine makeRandomPhase(field)
		implicit none
		type(fluid_field) :: field
		double precision :: r1, r2
		integer :: i, j, k, n
		integer, dimension(:), allocatable :: seed

		call RANDOM_SEED(size = n)
		allocate(seed(n))
		seed = RANK + RSEED
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
		seed = RANK + RSEED
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
		seed = RANK + RSEED
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

	subroutine insert_vortex_ring(field,x0,y0,z0,r0,dir,rot)
		implicit none
		integer :: i,j,k,dir,rot
		type(fluid_field) :: field
		double precision :: x0,y0,z0,r0
		double precision, dimension(sx:ex,sy:ey,sz:ez) :: rr1,rr2,d1,d2
		double precision, dimension(:,:), allocatable :: s

		if (rot .eq. 0) then
			allocate(s(sx:ex,sy:ey))
			!$OMP parallel do private (i,j) collapse(2)
			do j=sy,ey
				do i=sx,ex
					s(i,j) = sqrt((GX(i)-x0)**2 + (GY(j)-y0)**2)
				end do
			end do
			!$OMP end parallel do
		else if (rot .eq. 1) then
			allocate(s(sx:ex,sz:ez))
			!$OMP parallel do private (i,k) collapse(2)
			do k=sz,ez
				do i=sx,ex
					s(i,k) = sqrt((GX(i)-x0)**2 + (GZ(k)-z0)**2)
				end do
			end do
			!$OMP end parallel do
		else if (rot .eq. 2) then
			allocate(s(sy:ey,sz:ez))
			!$OMP parallel do private (k,j) collapse(2)
			do k=sz,ez
				do j=sy,ey
					s(j,k) = sqrt((GY(j)-y0)**2 + (GZ(k)-z0)**2)
				end do
			end do
			!$OMP end parallel do
		end if

		!$OMP parallel do private (i,j,k) collapse(3)
		do k=sz,ez
			do j=sy,ey
				do i=sx,ex
				if (rot .eq. 0) then
					d1(i,j,k) = sqrt((GZ(k) - z0)**2 + (s(i,j)+r0)**2)
					d2(i,j,k) = sqrt((GZ(k) - z0)**2 + (s(i,j)-r0)**2)
				else if (rot .eq. 1) then
					d1(i,j,k) = sqrt((GY(j) - y0)**2 + (s(i,k)+r0)**2)
					d2(i,j,k) = sqrt((GY(j) - y0)**2 + (s(i,k)-r0)**2)
				else if (rot .eq. 2) then
					d1(i,j,k) = sqrt((GX(i) - x0)**2 + (s(j,k)+r0)**2 )
					d2(i,j,k) = sqrt((GX(i) - x0)**2 + (s(j,k)-r0)**2 )
				end if
				end do
			end do
		end do
		!$OMP end parallel do

		rr1 = sqrt(((0.3437d0+0.0572d0*d1**2))/(1.0d0+(0.6666d0*d1**2)+(0.0572d0*d1**4)))
		rr2 = sqrt(((0.3437d0+0.0572d0*d2**2))/(1.0d0+(0.6666d0*d2**2)+(0.0572d0*d2**4)))

		!$OMP parallel do private (i,j,k) collapse(3)
		do k=sz,ez
			do j=sy,ey
				do i=sx,ex
					if (rot .eq. 0) then
					field%GRID(i,j,k) = field%GRID(i,j,k)*rr1(i,j,k)*(GZ(k)-z0+dir*EYE*(s(i,j)+r0))* &
														  rr2(i,j,k)*(GZ(k)-z0-dir*EYE*(s(i,j)-r0))
					else if (rot .eq. 1) then
					field%GRID(i,j,k) = field%GRID(i,j,k)*rr1(i,j,k)*(GY(j)-y0+dir*EYE*(s(i,k)+r0))* &
														  rr2(i,j,k)*(GY(j)-y0-dir*EYE*(s(i,k)-r0))
					else if (rot .eq. 2) then
					field%GRID(i,j,k) = field%GRID(i,j,k)*rr1(i,j,k)*(GX(i)-x0+dir*EYE*(s(j,k)+r0))* &
														  rr2(i,j,k)*(GX(i)-x0-dir*EYE*(s(j,k)-r0))
					end if
				end do
			end do
		end do
		!$OMP end parallel do
		deallocate(s)
	end subroutine
end module