module workspace
    use params
    use parallel
    complex*16, dimension(:,:,:), ALLOCATABLE, TARGET :: GRID_RK4,GRID_RK4_T1
	complex*16, dimension(:,:,:), ALLOCATABLE :: GRID_T2,GRID_T3,DDI_K
	double precision, dimension(:), ALLOCATABLE :: GX,GY,GZ,KX,KY,KZ
	integer :: sx,sy,sz,ex,ey,ez
	double precision, dimension(:,:,:), ALLOCATABLE :: POT
	complex(C_DOUBLE_COMPLEX), pointer :: GRID(:,:,:),GRID_T1(:,:,:)
	double precision :: TIME,DKSPACE,EDD_T1,EDD_T2
    contains
	subroutine init_workspace
		implicit none

        if(METHOD==0) then
        	call calc_local_idx_3DWithGhost(NX,NY,NZ,sx,ex,sy,ey,sz,ez)
			ALLOCATE(GRID_RK4(sx:ex,sy:ey,sz:ez))
			ALLOCATE(GRID_RK4_T1(sx:ex,sy:ey,sz:ez))
			ALLOCATE(GRID_T2(sx:ex,sy:ey,sz:ez))
			ALLOCATE(GRID_T3(sx:ex,sy:ey,sz:ez))
			GRID => GRID_RK4
			GRID_T1 => GRID_RK4_T1
		else if (METHOD==1) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Distributing the computational grid along z"
			end if
			call parallel_barrier
			call setup_local_allocation(NX,NY,NZ,GRID,GRID_T1)
			call parallel_barrier
			sx = 1
			ex = NX
			sy = 1
			ey = NY
			sz = 1
			ez = local_NZ
			EDD_T1 = 1.0d0/(1.0d0-EDD)
			EDD_T2 = EDD/(1.0d0-EDD)
        end if
		ALLOCATE(POT(sx:ex,sy:ey,sz:ez))
		ALLOCATE(GX(sx:ex))
		ALLOCATE(GY(sy:ey))
		ALLOCATE(GZ(sz:ez))
		ALLOCATE(KX(sx:ex))
		ALLOCATE(KY(sy:ey))
		ALLOCATE(KZ(sz:ez))
		if(RANK .eq. 0) then
			write(6, '(a)') "Generating grid coordinates."
		end if
		call parallel_barrier
        call setupGXYZ
        call setupKXYZ
        if (METHOD==1) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Calculating k-space dipole-dipole potential."
			end if
			call parallel_barrier
			ALLOCATE(DDI_K(sx:ex,sy:ey,sz:ez))
        	call setupDDIK
        end if
        GRID = 0.0d0
		TIME = 0.0d0
	end subroutine

	subroutine setupDDIK
	    implicit none
	    integer :: i, j ,k
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    DDI_K(i,j,k) = 3.0d0*(KZ(k)**2.0d0)/(KX(i)**2.0d0+KY(j)**2.0d0+KZ(k)**2.0d0+1e-15) - 1.0d0
                end do
            end do
        end do
        !$OMP end parallel do
	end subroutine

	subroutine setupGXYZ
	    implicit none
	    integer :: i, j ,k
		do k = sz,ez
		    GZ(k) = (k+local_k_offset-1)*DSPACE
		end do		

		do j = sy,ey
		    GY(j) = (j-1)*DSPACE
		end do

		do i = sx,ex
		    GX(i) = (i-1)*DSPACE
		end do
	end subroutine

	subroutine setupKXYZ
	    implicit none
	    integer :: i, j ,k
	    DKSPACE = PI/(DSPACE*NZ)

		do k = sz,ez
			if ((k+local_k_offset-1) <= NX/2) then
				KZ(k) = 2.0d0*PI*(k+local_k_offset-1)/(DSPACE*NZ)
			else
				KZ(k) = 2.0d0*PI*(NZ-k-local_k_offset+1)/(DSPACE*NZ)
			end if
		end do		

		do j = sy,ey
			if (j-1 <= NY/2) then
				KY(j) = 2.0d0*PI*(j-1)/(DSPACE*NY)
			else
				KY(j) = 2.0d0*PI*(NY-j+1)/(DSPACE*NY)
			end if
		end do

		do i = sx,ex
			if (i-1 <= NX/2) then
				KX(i) = 2.0d0*PI*(i-1)/(DSPACE*NX)
			else
				KX(i) = 2.0d0*PI*(NX-i+1)/(DSPACE*NX)
			end if
		end do

	end subroutine

end module