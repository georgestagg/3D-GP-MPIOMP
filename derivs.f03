module derivs
use workspace
interface laplacian
	procedure laplacian_fluid, laplacian_magnetic
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
	double precision function laplacian_magnetic(field_in,i,j,k)
		use params
		implicit none
		class(magnetic_field) :: field_in
		integer :: i,j,k
		laplacian_magnetic = (-6.0d0*field_in%GRID_R(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k)&
								  + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k)&
								  + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
	end function

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
			d2dx2 = (-2.0d0*field_in%GRID_R(i,j,k) + BC(field_in,i+1,j,k)+BC(field_in,i-1,j,k))/(DSPACE**2.0d0)
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
			d2dy2 = (-2.0d0*field_in%GRID_R(i,j,k) + BC(field_in,i,j+1,k)+BC(field_in,i,j-1,k))/(DSPACE**2.0d0)
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
			d2dz2 = (-2.0d0*field_in%GRID_R(i,j,k) + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
		class is (fluid_field)
			d2dz2 = (-2.0d0*field_in%GRID(i,j,k) + BC(field_in,i,j,k+1)+BC(field_in,i,j,k-1))/(DSPACE**2.0d0)
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
			BC = field_in%GRID_R(ii,jj,kk)
		class is (fluid_field)
			BC = field_in%GRID(ii,jj,kk)
			!Quasi-periodic fluid field
			if(i == NX+1 .and. BCX==3) then
			  BC = BC*exp(EYE*PI*NVORTALL(3,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
			              -EYE*PI*NVORTALL(2,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(i == 0 .and. BCX==3) then
			  BC =BC/exp(EYE*PI*NVORTALL(3,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
			             -EYE*PI*NVORTALL(2,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(j == NY+1 .and. BCY==3) then
			  BC = BC*exp(-EYE*PI*NVORTALL(3,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
			              +EYE*PI*NVORTALL(1,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(j == 0 .and. BCY==3) then
			  BC = BC/exp(-EYE*PI*NVORTALL(3,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
			              +EYE*PI*NVORTALL(1,field_in%field_number)*GZ(kk)/((NZ-1)*DSPACE))
			end if
			if(k == NZ+1 .and. BCZ==3) then
			  BC = BC*exp(EYE*PI*NVORTALL(2,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
			              -EYE*PI*NVORTALL(1,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
						  +EYE*PI*NVORTALL(1,field_in%field_number)*NVORTALL(2,field_in%field_number))
			end if
			if(k == 0 .and. BCZ==3) then
			  BC = BC/exp(EYE*PI*NVORTALL(2,field_in%field_number)*GX(ii)/((NX-1)*DSPACE)&
			              -EYE*PI*NVORTALL(1,field_in%field_number)*GY(jj)/((NY-1)*DSPACE)&
						  +EYE*PI*NVORTALL(1,field_in%field_number)*NVORTALL(2,field_in%field_number))
			end if
		end select
	end function
end module