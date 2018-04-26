module workspace
    use params
    use parallel

	type :: computational_field
		complex(C_DOUBLE_COMPLEX), pointer :: GRID(:,:,:)
		complex(C_DOUBLE_COMPLEX), pointer :: GRID_T(:,:,:)
	end type computational_field
	type(computational_field),DIMENSION(:),POINTER :: WS

	complex(C_DOUBLE_COMPLEX), pointer :: GRID_T1(:,:,:),GRID_T2(:,:,:)
	complex*16, dimension(:,:,:), ALLOCATABLE :: DDI_K,GRID_E,QP_EX,QP_EY,QP_EZ
	double precision, dimension(:), ALLOCATABLE :: GX,GY,GZ,KX,KY,KZ
	double precision, dimension(:,:), ALLOCATABLE :: GG
	integer :: sx,sy,sz,ex,ey,ez,cur_step
	double precision, dimension(:,:,:), ALLOCATABLE :: POT
	double precision :: TIME,DKSPACE,EDD_T1,EDD_T2
    contains
	subroutine init_workspace
		implicit none
	    integer :: f
		if(RANK .eq. 0) then
			write(6, '(a)') "Allocating distributed workspace data"
		end if

		ALLOCATE(WS(FIELDS))
		ALLOCATE(GG(FIELDS,FIELDS))

        if(METHOD==0) then
        	call calc_local_idx_3DWithGhost(NX,NY,NZ,sx,ex,sy,ey,sz,ez)
			do f = 1,FIELDS
				ALLOCATE(WS(f)%GRID(sx:ex,sy:ey,sz:ez))
				ALLOCATE(WS(f)%GRID_T(sx:ex,sy:ey,sz:ez))
				WS(f)%GRID = 0.0d0
				WS(f)%GRID_T = 0.0d0
			end do
			ALLOCATE(GRID_T1(sx:ex,sy:ey,sz:ez))
			ALLOCATE(GRID_T2(sx:ex,sy:ey,sz:ez))
		else if (METHOD==1) then
			call setup_local_allocation(NX,NY,NZ,WS(1)%GRID,GRID_T1)
			call parallel_barrier
			sx = 1
			ex = NX
			sy = 1
			ey = NY
			sz = 1
			ez = local_NZ
			EDD_T1 = 1.0d0/(1.0d0-EDD)
			EDD_T2 = EDD/(1.0d0-EDD)
			WS(1)%GRID = 0.0d0
        end if
		ALLOCATE(POT(sx:ex,sy:ey,sz:ez))
		ALLOCATE(GX(sx:ex))
		ALLOCATE(GY(sy:ey))
		ALLOCATE(GZ(sz:ez))
		ALLOCATE(KX(sx:ex))
		ALLOCATE(KY(sy:ey))
		ALLOCATE(KZ(sz:ez))
		if(RANK .eq. 0) then
			write(6, '(a)') "Calculating position and momentum space coordinates"
		end if
		call parallel_barrier
        call setupGXYZ
        call setupKXYZ
        if (METHOD==1) then
			if(RANK .eq. 0) then
				write(6, '(a)') "Pre-calculating dipole-dipole k-space potential"
				write(6, '(a,f6.3)') "EDD is: ", EDD
			end if
			call parallel_barrier
			ALLOCATE(DDI_K(sx:ex,sy:ey,sz:ez))
        	call setupDDIK
        end if
		if (RHSType .eq. 3) then
			ALLOCATE(GRID_E(sx:ex,sy:ey,sz:ez))
			ALLOCATE(QP_EX(sx:ex,sy:ey,sz:ez))
			ALLOCATE(QP_EY(sx:ex,sy:ey,sz:ez))
			ALLOCATE(QP_EZ(sx:ex,sy:ey,sz:ez))
			call setupEXEYEZ
		end if
        GRID_T1 = 0.0d0
		TIME = 0.0d0
	end subroutine

	subroutine setupDDIK
	    implicit none
	    integer :: i, j ,k
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    !NOTE: This term is calculated in transposed k-space, so KY is really KZ.
                    !NB: http://fftw.org/doc/Transposed-distributions.html#Transposed-distributions
                    DDI_K(i,j,k) = 3.0d0*(KY(j)**2.0d0)/(KX(i)**2.0d0+KY(j)**2.0d0+KZ(k)**2.0d0+1e-20) - 1.0d0
                end do
            end do
        end do
        !$OMP end parallel do
	end subroutine

	subroutine setupEXEYEZ
	    implicit none
	    integer :: i, j ,k
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    QP_EX(i,j,k) = exp(EYE*PI*NVORTY*GX(i)*GZ(k)/((NX-1)*DSPACE*(NZ-1)*DSPACE)&
                                  -EYE*PI*NVORTZ*GX(i)*GY(j)/((NX-1)*DSPACE*(NY-1)*DSPACE))
                end do
            end do
        end do
        !$OMP end parallel do
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    QP_EY(i,j,k) = exp(-EYE*PI*NVORTX*GY(j)*GZ(k)/((NY-1)*DSPACE*(NZ-1)*DSPACE)&
                                  +EYE*PI*NVORTZ*GY(j)*GX(i)/((NY-1)*DSPACE*(NX-1)*DSPACE))
                end do
            end do
        end do
        !$OMP end parallel do
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    QP_EZ(i,j,k) = exp(EYE*PI*NVORTX*GZ(k)*GY(j)/((NZ-1)*DSPACE*(NY-1)*DSPACE)&
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

end module