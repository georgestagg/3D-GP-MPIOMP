module io
	use workspace
	use netcdf
	implicit none
	include 'mpif.h'
	integer, allocatable :: f_re_id(:), f_im_id(:)
	integer :: ncdf_id,x_id,y_id,z_id,pot_id,time_id,step_id,r,mpi_info, stat
	integer       :: icount(3)
	integer       :: istarting(3)
	Contains

	subroutine dump_wavefunction (II)
		implicit none
		integer :: II
		character(len=30) :: fname
		if(RT == 1) then
			write(fname, '(a,i0.6,a)') 'psi.', II/DUMPWF,".nc"
		else
			write(fname, '(a,i0.6,a)') 'imag.', II/DUMPWF,".nc"
		end if
		if(RANK .eq. 0) then
			write(6,'(a,a)',advance="no") "Writing: ", TRIM(fname)
		end if

		call make_file(TRIM(fname))

		if(RANK .eq. 0) then
			write(6,'(a)',advance="no") "."
		end if

		call write_wf_file(rt)

		if(RANK .eq. 0) then
			write(6,'(a)',advance="no") "."
		end if

		call close_file

		if(RANK .eq. 0) then
			write(6,'(a)',advance="no") "."
			write(6,'(a)') "Done!"
		end if
	end subroutine

	subroutine make_file(fname)
		implicit none
		integer :: dims(3), f
		integer :: x_dim_id,y_dim_id,z_dim_id
		character(len=*) :: fname
		character(len=1024) :: vname

		call MPI_Info_create(mpi_info, IERR)
		r=NF90_create(fname, IOR(NF90_NETCDF4, NF90_MPIIO), ncdf_id)
		call handle_err(r)

		r=NF90_def_dim(ncdf_id, 'x_dim', NX, x_dim_id)
		call handle_err(r)
		r=NF90_def_dim(ncdf_id, 'y_dim', NY, y_dim_id)
		call handle_err(r)
		r=NF90_def_dim(ncdf_id, 'z_dim', NZ, z_dim_id)
		call handle_err(r)

		r=NF90_def_var(ncdf_id, 'gx', NF90_DOUBLE,x_dim_id, x_id)
		call handle_err(r)
		r=NF90_def_var(ncdf_id, 'gy', NF90_DOUBLE,y_dim_id, y_id)
		call handle_err(r)
		r=NF90_def_var(ncdf_id, 'gz', NF90_DOUBLE,z_dim_id, z_id)
		call handle_err(r)

		dims(1) = x_dim_id
		dims(2) = y_dim_id
		dims(3) = z_dim_id

		deallocate(f_re_id,STAT=stat)
		allocate(f_re_id(FLUIDS))
		deallocate(f_im_id,STAT=stat)
		allocate(f_im_id(FLUIDS))

		do f = 1,FLUIDS
			write(vname,'(a,i0.3,a)') 'fluid_',f,'_real'
			r=NF90_def_var(ncdf_id, vname, NF90_DOUBLE, dims, f_re_id(f))
			call handle_err(r)
			write(vname,'(a,i0.3,a)') 'fluid_',f,'_imag'
			r=NF90_def_var(ncdf_id, vname, NF90_DOUBLE, dims, f_im_id(f))
			call handle_err(r)
		end do

		r=NF90_def_var(ncdf_id, 'pot' , NF90_DOUBLE, dims, pot_id)
		call handle_err(r)

		r=NF90_def_var(ncdf_id, 'step' , NF90_INT, step_id)
		call handle_err(r)
		r=NF90_def_var(ncdf_id, 'time' , NF90_DOUBLE, time_id)
		call handle_err(r)

		r = NF90_enddef(ncdf_id)
		call handle_err(r)
	end subroutine

	subroutine write_wf_file(rt)
		implicit none
		integer :: rt
		if(METHOD==0) then
			call write_wf_file_RK4(rt)
		else if(METHOD==1) then
			call write_wf_file_FFTW(rt)
		end if
	end subroutine

	subroutine write_wf_file_RK4(rt)
		implicit none
		integer :: rt,f

		r = nf90_inq_varid(ncdf_id, "gx", x_id)
		call handle_err(r)
		!Dont forget to ignore ghost points!
		istarting(1) = sx+NGHOST
		icount(1) = ex-sx-2*NGHOST+1
		r=NF90_put_var(ncdf_id, x_id, GX(sx+NGHOST:ex-NGHOST),istarting,icount)
		call handle_err(r)
		istarting(1) = sy+NGHOST
		icount(1) = ey-sy-2*NGHOST+1
		r=NF90_put_var(ncdf_id, y_id, GY(sy+NGHOST:ey-NGHOST),istarting,icount)
		call handle_err(r)
		istarting(1) = sz+NGHOST
		icount(1) = ez-sz-2*NGHOST+1
		r=NF90_put_var(ncdf_id, z_id, GZ(sz+NGHOST:ez-NGHOST),istarting,icount)
		call handle_err(r)

		istarting(1) = sx+NGHOST
		icount(1) = ex-sx-2*NGHOST+1
		istarting(2) = sy+NGHOST
		icount(2) = ey-sy-2*NGHOST+1
		istarting(3) = sz+NGHOST
		icount(3) = ez-sz-2*NGHOST+1
		do f = 1,FLUIDS
			r=NF90_put_var(ncdf_id,f_re_id(f),DBLE(WS%FLUID(f)%GRID(sx+NGHOST:ex-NGHOST, &
						sy+NGHOST:ey-NGHOST,sz+NGHOST:ez-NGHOST)),istarting,icount)
			call handle_err(r)
			r=NF90_put_var(ncdf_id,f_im_id(f),DIMAG(WS%FLUID(f)%GRID(sx+NGHOST:ex-NGHOST, &
						sy+NGHOST:ey-NGHOST,sz+NGHOST:ez-NGHOST)),istarting,icount)
			call handle_err(r)
		end do

		r=NF90_put_var(ncdf_id,pot_id,POT(sx+NGHOST:ex-NGHOST,sy+NGHOST:ey-NGHOST,sz+NGHOST:ez-NGHOST),istarting,icount)
		call handle_err(r)

		if(rt == 1) then
			r=NF90_put_var(ncdf_id,step_id,cur_step)
		else
			r=NF90_put_var(ncdf_id,step_id,0)
		end if
		call handle_err(r)

		r=NF90_put_var(ncdf_id,time_id,TIME)
		call handle_err(r)
	end subroutine

	subroutine write_wf_file_FFTW(rt)
		implicit none
		integer :: rt
		istarting(1) = 1
		icount(1) = NX
		r=NF90_put_var(ncdf_id,x_id,GX,istarting,icount)
		call handle_err(r)
		istarting(1) = 1
		icount(1) = NY
		r=NF90_put_var(ncdf_id,y_id,GY,istarting,icount)
		call handle_err(r)
		istarting(1) = 1+local_k_offset
		icount(1) = ez
		r=NF90_put_var(ncdf_id,z_id,GZ,istarting,icount)
		call handle_err(r)

		istarting(1) = 1
		icount(1) = NX
		istarting(2) = 1
		icount(2) = NY
		istarting(3) = 1+local_k_offset
		icount(3) = ez
		r=NF90_put_var(ncdf_id,f_re_id(1),DBLE(WS%FLUID(1)%GRID),istarting,icount)
		call handle_err(r)
		r=NF90_put_var(ncdf_id,f_im_id(1),DIMAG(WS%FLUID(1)%GRID),istarting,icount)
		call handle_err(r)

		r=NF90_put_var(ncdf_id,pot_id,POT,istarting,icount)
		call handle_err(r)

		if(rt == 1) then
			r=NF90_put_var(ncdf_id,step_id,cur_step)
		else
			r=NF90_put_var(ncdf_id,step_id,0)
		end if
		call handle_err(r)
		r=NF90_put_var(ncdf_id,time_id,TIME)
		call handle_err(r)
	end subroutine

	subroutine close_file()
		implicit none
		r=NF90_sync(ncdf_id)
		call handle_err(r)
		r=NF90_close(ncdf_id)
		call handle_err(r)
	end subroutine

	subroutine read_wf_file(fname)
		implicit none
		character(len=2048) fname
		if(METHOD==0) then
			call read_wf_file_RK4(fname)
		else if(METHOD==1) then
			call read_wf_file_FFTW(fname)
		end if
	end subroutine

	subroutine read_wf_file_RK4(fname)
		implicit none
		integer :: rwf_ncid,rwf_pot_id,rwf_step_id,rwf_time_id,f
		integer, allocatable :: rwf_re_id(:),rwf_im_id(:)
		double precision, dimension(sx:ex,sy:ey,sz:ez) :: realgrid,imaggrid
		character(len=*) fname
		character(len=1024) vname
		if(RANK .eq. 0) then
			write(6,*) "Reading saved state..."
			write(6,*) "Opening file: ", TRIM(fname)
		end if

		r = NF90_open(fname,NF90_NOWRITE,rwf_ncid)
		call handle_err(r)

		deallocate(rwf_re_id,STAT=stat)
		allocate(rwf_re_id(FLUIDS))
		deallocate(rwf_im_id,STAT=stat)
		allocate(rwf_im_id(FLUIDS))

		do f = 1,FLUIDS
			write(vname,'(a,i0.3,a)') 'fluid_',f,'_real'
			r = NF90_inq_varid(rwf_ncid, vname,  rwf_re_id(f))
			call handle_err(r)
			write(vname,'(a,i0.3,a)') 'fluid_',f,'_imag'
			r = NF90_inq_varid(rwf_ncid, vname,  rwf_im_id(f))
			call handle_err(r)
		end do

		r = NF90_inq_varid(rwf_ncid,  "pot", rwf_pot_id)
		call handle_err(r)
		r = NF90_inq_varid(rwf_ncid, "step",  rwf_step_id)
		call handle_err(r)
		r = NF90_inq_varid(rwf_ncid,  "time", rwf_time_id)
		call handle_err(r)

		istarting(1) = sx+1
		icount(1) = ex-sx-1
		istarting(2) = sy+1
		icount(2) = ey-sy-1
		istarting(3) = sz+1
		icount(3) = ez-sz-1
		do f = 1,FLUIDS
			r = NF90_get_var(rwf_ncid, rwf_re_id(f), realgrid(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),start=istarting,count=icount)
			call handle_err(r)
			r = NF90_get_var(rwf_ncid, rwf_im_id(f), imaggrid(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),start=istarting,count=icount)
			call handle_err(r)
			WS%FLUID(f)%GRID = realgrid + EYE*imaggrid
		end do

		r = NF90_get_var(rwf_ncid, rwf_pot_id, POT(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),start=istarting,count=icount)
		call handle_err(r)

		!r=NF90_get_var(rwf_ncid,rwf_step_id,INITSSTEP)
		!call handle_err(r)
		r=NF90_get_var(rwf_ncid,rwf_time_id,TIME)
		call handle_err(r)

		r = NF90_close(rwf_ncid)
		call handle_err(r)

		if(RANK .eq. 0) then
			write(6,*) "Reload complete!"
			write(6,fmt="(a,i8,a,f10.3)") "Resuming simulation from step number: ",INITSTEP, ", at time: ", TIME
		end if
	end subroutine


	subroutine read_surf_file_RK4(fname)
		implicit none
		integer :: rsurf_ncid,rsurf_surf_id,i,j,k
		double precision, dimension(sx:ex,sy:ey) :: h
		double precision :: localhmin, globalhmin
		character(len=2048) fname
		if(RANK .eq. 0) then
			write(6,*) "Opening surface file: ", TRIM(fname)
		end if

		h = 0.0d0
		
		r = NF90_open(fname,IOR(nf90_netcdf4,nf90_MPIIO),rsurf_ncid)
		call handle_err(r)
		r = NF90_inq_varid(rsurf_ncid,  "surf", rsurf_surf_id)
		call handle_err(r)

		istarting(1) = sx+1
		icount(1) = ex-sx-1
		istarting(2) = sy+1
		icount(2) = ey-sy-1

		r = NF90_get_var(rsurf_ncid, rsurf_surf_id, h(sx+1:ex-1,sy+1:ey-1),start=istarting,count=icount)
		call handle_err(r)

		r = NF90_close(rsurf_ncid)
		call handle_err(r)

		if(RANK .eq. 0) then
			write(6,*) "Surface file read complete!"
		end if 

		localhmin = minval(h)
		call MPI_Allreduce(localhmin, globalhmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_GRID,IERR_3DWG)

		!$OMP parallel do private (i,j,k) collapse(3)
	    do k = sz+1,ez-1
	        do j = sy+1,ey-1
	            do i = sx+1,ex-1
	            	if (GZ(k) < h(i,j) - globalhmin) then
	                	POT(i,j,k) = OBJHEIGHT
	                else
	                	POT(i,j,k) = 0.0d0
	                end if
	            end do
	        end do
	    end do
    	!$OMP end parallel do

		if(RANK .eq. 0) then
			write(6,"(a)") "Potential successfully built from surface."
		end if
	end subroutine


	subroutine read_wf_file_FFTW(fname)
		implicit none
		integer :: rwf_ncid,rwf_pot_id,rwf_step_id,rwf_time_id
		double precision, dimension(sx:ex,sy:ey,sz:ez) :: realgrid,imaggrid
		integer, allocatable :: rwf_re_id(:),rwf_im_id(:)
		character(len=*) fname
		if(RANK .eq. 0) then
			write(6,*) "Reading saved state..."
			write(6,*) "Opening file: ", TRIM(fname)
		end if

		r = NF90_open(fname,IOR(nf90_netcdf4,nf90_MPIIO),rwf_ncid)
		call handle_err(r)

		deallocate(rwf_re_id,STAT=stat)
		allocate(rwf_re_id(1))
		deallocate(rwf_im_id,STAT=stat)
		allocate(rwf_im_id(1))

		r = NF90_inq_varid(rwf_ncid, 'fluid_001_real',  rwf_re_id(1))
		call handle_err(r)
		r = NF90_inq_varid(rwf_ncid, 'fluid_001_imag',  rwf_im_id(1))
		call handle_err(r)

		r = NF90_inq_varid(rwf_ncid,  "pot", rwf_pot_id)
		call handle_err(r)
		r = NF90_inq_varid(rwf_ncid, "step",  rwf_step_id)
		call handle_err(r)
		r = NF90_inq_varid(rwf_ncid, "time", rwf_time_id)
		call handle_err(r)

		istarting(1) = 1
		icount(1) = NX
		istarting(2) = 1
		icount(2) = NY
		istarting(3) = 1+local_k_offset
		icount(3) = ez

		r = NF90_get_var(rwf_ncid, rwf_re_id(1), realgrid,start=istarting,count=icount)
		call handle_err(r)
		r = NF90_get_var(rwf_ncid, rwf_im_id(1), imaggrid,start=istarting,count=icount)
		call handle_err(r)
		WS%FLUID(1)%GRID = realgrid + EYE*imaggrid

		r = NF90_get_var(rwf_ncid, rwf_pot_id, POT,start=istarting,count=icount)
		call handle_err(r)

		!r=NF90_get_var(rwf_ncid,rwf_step_id,INITSSTEP)
		!call handle_err(r)
		r=NF90_get_var(rwf_ncid,rwf_time_id,TIME)
		call handle_err(r)

		r = NF90_close(rwf_ncid)
		call handle_err(r)

		if(RANK .eq. 0) then
			write(6,*) "Reload complete!"
			write(6,fmt="(a,i8,a,f10.3)") "Resuming simulation from step number: ",INITSTEP, ", at time: ", TIME
		end if
	end subroutine

	subroutine handle_err(status)
		integer, intent ( in) :: status
		if(status /= nf90_noerr) then
			print *, r, trim(nf90_strerror(status))
			call abort()
		end if
	end subroutine handle_err

end module
