module workspace
    use params
    use parallel
    complex*16, dimension(:,:,:), ALLOCATABLE, TARGET :: GRID_RK4
	complex*16, dimension(:,:,:), ALLOCATABLE :: GRID_T1,GRID_T2,GRID_T3
	double precision, dimension(:), ALLOCATABLE :: GX,GY,GZ
	integer :: sx,sy,sz,ex,ey,ez
	double precision, dimension(:,:,:), ALLOCATABLE :: POT
	complex(C_DOUBLE_COMPLEX), pointer :: GRID(:,:,:)
	double precision :: TIME
    contains
	subroutine init_workspace
		implicit none

        if(METHOD==0) then
        	call calc_local_idx_3DWithGhost(NX,NY,NZ,sx,ex,sy,ey,sz,ez)
            ALLOCATE(GRID_RK4(sx:ex,sy:ey,sz:ez))
			GRID => GRID_RK4
		else if (METHOD==1) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Distributing the computational grid along z"
			end if
			call parallel_barrier
			call setup_local_allocation(NX,NY,NZ,GRID)
			call parallel_barrier
			sx = 1
			ex = NX
			sy = 1
			ey = NY
			sz = 1
			ez = local_NZ
        end if
		ALLOCATE(POT(sx:ex,sy:ey,sz:ez))
		ALLOCATE(GX(sx:ex))
		ALLOCATE(GY(sy:ey))
		ALLOCATE(GZ(sz:ez))
		ALLOCATE(GRID_T1(sx:ex,sy:ey,sz:ez))
		ALLOCATE(GRID_T2(sx:ex,sy:ey,sz:ez))
		ALLOCATE(GRID_T3(sx:ex,sy:ey,sz:ez))
		if(RANK .eq. 0) then
			write(6, '(a)') "Generating grid coordinates."
		end if
		call parallel_barrier
        call setupGXYZ
        GRID = 0.0d0
		TIME = 0.0d0
	end subroutine

	subroutine setupGXYZ
	    implicit none
	    integer :: i, j ,k
		do i = sx,ex
		    GX(i) = (i-1)*DSPACE
		end do
		
		do j = sy,ey
		    GY(j) = (j-1)*DSPACE
		end do
		
		do k = sz,ez
		    GZ(k) = (k+local_k_offset-1)*DSPACE
		end do
	end subroutine

end module