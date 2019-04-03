module derivs
use workspace
interface laplacian
	procedure laplacian_fluid
end interface
interface BC
	procedure MBC,FBC
end interface
contains
	complex*16 function laplacian_fluid(field_in,i,j,k)
		use params
		implicit none
		class(fluid_field) :: field_in
		integer :: i,j,k
		laplacian_fluid = (-6.0d0*field_in%GRID(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k)&
								  + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k)&
								  + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
	end function

	COMPLEX*16 function ddx(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(fluid_field) :: field_in
		ddx = (BC(field_in,i,j,k)-BC(field_in,i-1,j,k))/DSPACE
	end function


	COMPLEX*16 function ddy(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(fluid_field) :: field_in
		ddy = (BC(field_in,i,j,k)-BC(field_in,i,j-1,k))/DSPACE
	end function

	COMPLEX*16 function ddz(field_in,i,j,k)
		use params
		implicit none
		class(fluid_field) :: field_in
		integer :: i,j,k
		ddz = (BC(field_in,i,j,k)-BC(field_in,i,j,k-1))/DSPACE
	end function

	COMPLEX*16 function d2dx2(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(fluid_field) :: field_in
		d2dx2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k))/(DSPACE**2.0d0)
	end function

	COMPLEX*16 function d2dy2(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(fluid_field) :: field_in
		d2dy2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k))/(DSPACE**2.0d0)
	end function

	COMPLEX*16 function d2dz2(field_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		class(fluid_field) :: field_in
		d2dz2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
	end function

	double precision function curlcurlx_magnetic(ws_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		type(workspace_t), intent(in) :: ws_in

		curlcurlx_magnetic = (4.0d0*BC(ws_in,i,j,k,1) - BC(ws_in,i,j+1,k,1) - BC(ws_in,i,j-1,k,1) - BC(ws_in,i,j,k+1,1) &
							- BC(ws_in,i,j,k-1,1) &
							+ BC(ws_in,i+1,j,k,2) - BC(ws_in,i,j,k,2) - BC(ws_in,i+1,j-1,k,2) + BC(ws_in,i,j-1,k,2) &
							+ BC(ws_in,i+1,j,k,3) - BC(ws_in,i,j,k,3) - BC(ws_in,i+1,j,k-1,3) + BC(ws_in,i,j,k-1,3))/(DSPACE**2.0d0)
	end function

	double precision function curlcurly_magnetic(ws_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		type(workspace_t), intent(in) :: ws_in

		curlcurly_magnetic = (4.0d0*BC(ws_in,i,j,k,2) - BC(ws_in,i,j,k+1,2) - BC(ws_in,i,j,k-1,2) - BC(ws_in,i+1,j,k,2) &
							- BC(ws_in,i-1,j,k,2) &
							+ BC(ws_in,i,j+1,k,3) - BC(ws_in,i,j,k,3) - BC(ws_in,i,j+1,k-1,3) + BC(ws_in,i,j,k-1,3) &
							+ BC(ws_in,i,j+1,k,1) - BC(ws_in,i,j,k,1) - BC(ws_in,i-1,j+1,k,1) + BC(ws_in,i-1,j,k,1))/(DSPACE**2.0d0)
	end function

	double precision function curlcurlz_magnetic(ws_in,i,j,k)
		use params
		implicit none
		integer :: i,j,k
		type(workspace_t), intent(in) :: ws_in

		curlcurlz_magnetic = (4.0d0*BC(ws_in,i,j,k,3) - BC(ws_in,i+1,j,k,3) - BC(ws_in,i-1,j,k,3) - BC(ws_in,i,j+1,k,3) &
							- BC(ws_in,i,j-1,k,3) &
							+ BC(ws_in,i,j,k+1,1) - BC(ws_in,i,j,k,1) - BC(ws_in,i-1,j,k+1,1) + BC(ws_in,i-1,j,k,1) &
							+ BC(ws_in,i,j,k+1,2) - BC(ws_in,i,j,k,2) - BC(ws_in,i,j-1,k+1,2) + BC(ws_in,i,j-1,k,2))/(DSPACE**2.0d0)
	end function

	double precision function MBC(ws_in,i,j,k,dir)
		use params
		implicit none
		integer :: dir,i,j,k
		type(workspace_t), intent(in) :: ws_in

!		if((i==0 .or. j==0 .or. k==0) .and. BCMX==0) then
!			select case (dir)
!			case(1)
!			MBC=H(3)*DSPACE - ws_in%MAGNETIC(2)%GRID_R(i+1,j,k) + ws_in%MAGNETIC(2)%GRID_R(i,j,k) + ws_in%MAGNETIC(1)%GRID_R(i,j+1,k)
!			case(2)
!			MBC=H(1)*DSPACE - ws_in%MAGNETIC(3)%GRID_R(i,j+1,k) + ws_in%MAGNETIC(3)%GRID_R(i,j,k) + ws_in%MAGNETIC(2)%GRID_R(i,j,k+1)
!			case(3)
!			MBC=H(2)*DSPACE - ws_in%MAGNETIC(1)%GRID_R(i,j,k+1) + ws_in%MAGNETIC(1)%GRID_R(i,j,k) + ws_in%MAGNETIC(3)%GRID_R(i+1,j,k)
!			end select
!			RETURN
!		end if
!
!		if((i==NX+1 .or. j==NY+1 .or. k==NZ+1) .and. BCMX==0) then
!			select case (dir)
!			case(1)
!			MBC= -H(3)*DSPACE + ws_in%MAGNETIC(2)%GRID_R(i,j,k) - ws_in%MAGNETIC(2)%GRID_R(i-1,j,k) + ws_in%MAGNETIC(1)%GRID_R(i,j-1,k)
!			case(2)
!			MBC= -H(1)*DSPACE + ws_in%MAGNETIC(3)%GRID_R(i,j,k) - ws_in%MAGNETIC(3)%GRID_R(i,j-1,k) + ws_in%MAGNETIC(2)%GRID_R(i,j,k-1)
!			case(3)
!			MBC= -H(2)*DSPACE + ws_in%MAGNETIC(1)%GRID_R(i,j,k) - ws_in%MAGNETIC(1)%GRID_R(i,j,k-1) + ws_in%MAGNETIC(3)%GRID_R(i-1,j,k)
!			end select
!			RETURN
!		end if

		MBC = ws_in%MAGNETIC(dir)%GRID_R(i,j,k)

	end function

	complex*16 function FBC(field_in,i,j,k,dir)
		use params
		implicit none
		class(fluid_field) :: field_in
		integer,optional :: dir
		integer :: i,j,k,ii,jj,kk
		!Zero
		if((i>=NX .or. i<=1) .and. BCX==2) then
			FBC = 0.0d0
			RETURN
		else if ((j>=NY .or. j<=1) .and. BCY==2) then
			FBC = 0.0d0
			RETURN
		else if ((k>=NZ .or. k<=1) .and. BCZ==2) then
			FBC = 0.0d0
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

		FBC = field_in%GRID(ii,jj,kk)
		!Quasi-periodic fluid field
		if(i == NX+1 .and. BCX==3) then
		  FBC = FBC*exp(EYE*PI*NVORTALL(3,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
		              -EYE*PI*NVORTALL(2,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
		end if
		if(i == 0 .and. BCX==3) then
		  FBC =FBC/exp(EYE*PI*NVORTALL(3,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
		             -EYE*PI*NVORTALL(2,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
		end if
		if(j == NY+1 .and. BCY==3) then
		  FBC = FBC*exp(-EYE*PI*NVORTALL(3,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
		              +EYE*PI*NVORTALL(1,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
		end if
		if(j == 0 .and. BCY==3) then
		  FBC = FBC/exp(-EYE*PI*NVORTALL(3,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
		              +EYE*PI*NVORTALL(1,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
		end if
		if(k == NZ+1 .and. BCZ==3) then
		  FBC = FBC*exp(EYE*PI*NVORTALL(2,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
		              -EYE*PI*NVORTALL(1,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
					  +EYE*PI*NVORTALL(1,field_in%field_number)*NVORTALL(2,field_in%field_number))
		end if
		if(k == 0 .and. BCZ==3) then
		  FBC = FBC/exp(EYE*PI*NVORTALL(2,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
		              -EYE*PI*NVORTALL(1,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
					  +EYE*PI*NVORTALL(1,field_in%field_number)*NVORTALL(2,field_in%field_number))
		end if
	end function
end module
