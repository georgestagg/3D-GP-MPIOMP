module output
    use workspace
    use netcdf
    implicit none
    include 'mpif.h'
    integer :: ncdf_id,x_id,y_id,z_id,re_id,im_id,pot_id,info
    integer       :: icount(3)
    integer       :: istarting(3)
    Contains

    subroutine dump_wavefunction (II,rt)
        implicit none
        integer :: II,rt
        character(len=80) fname
        if(rt .eq. 1) then
            write(fname, '(a,i0.6,a)') 'psi.', II/DUMPWF,".nc"
        else
            write(fname, '(a,i0.6,a)') 'imag.', II/DUMPWF,".nc"
        end if
        if(RANK .eq. 0) then
            write(6,'(a,a)',advance="no") "Writing: ", TRIM(fname)
        end if

        call make_file(fname)

        if(RANK .eq. 0) then
            write(6,'(a)',advance="no") "."
        end if

        call write_wf_file

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
        integer :: dims(3)
        integer :: r,x_dim_id,y_dim_id,z_dim_id
        character(len=*) :: fname

        call MPI_Info_create(info, IERR)
        call MPI_Info_set(info,"IBM_largeblock_io","true", IERR)

        r=NF90_create_par(fname , IOR(nf90_netcdf4, nf90_MPIIO), MPI_COMM, MPI_INFO_NULL,ncdf_id)
        
        call handle_err(r)
        r=NF90_def_dim(ncdf_id, 'x_dim', NX, x_dim_id)
        r=NF90_def_dim(ncdf_id, 'y_dim', NY, y_dim_id)
        r=NF90_def_dim(ncdf_id, 'z_dim', NZ, z_dim_id)

        r=NF90_def_var(ncdf_id, 'gx', NF90_DOUBLE,x_dim_id, x_id)
        r=NF90_def_var(ncdf_id, 'gy', NF90_DOUBLE,y_dim_id, y_id)
        r=NF90_def_var(ncdf_id, 'gz', NF90_DOUBLE,z_dim_id, z_id)

        dims(1) = x_dim_id
        dims(2) = y_dim_id
        dims(3) = z_dim_id

        r=NF90_def_var(ncdf_id, 'real', NF90_DOUBLE, dims, re_id)
        r=NF90_def_var(ncdf_id, 'imag', NF90_DOUBLE, dims, im_id)
        r=NF90_def_var(ncdf_id, 'pot' , NF90_DOUBLE, dims, pot_id)

        r=nf90_var_par_access(ncdf_id, x_id, NF90_COLLECTIVE)
        call handle_err(r)
    end subroutine

    subroutine write_wf_file
        implicit none
        if(METHOD==0) then
            call write_wf_file_RK4
        else if(METHOD==1) then
            call write_wf_file_FFTW
        end if
    end subroutine

    subroutine write_wf_file_RK4
        implicit none
        integer :: r,i
        !Dont forget to ignore ghost points!
        istarting(1) = sx+1
        icount(1) = ex-sx-1
        !write(6,'(a,i4,a,i3,a,i3,a,i3,a,i3)') "Rank: ", rank, " is writing GX(", sx+1,":",ex-1, &
        !    ") starting from",istarting(1),"with length",icount(1)
        r=NF90_put_var(ncdf_id, x_id, GX(sx+1:ex-1),istarting,icount)
        istarting(1) = sy+1
        icount(1) = ey-sy-1
        r=NF90_put_var(ncdf_id, y_id, GY(sy+1:ey-1),istarting,icount)
        istarting(1) = sz+1
        icount(1) = ez-sz-1
        r=NF90_put_var(ncdf_id, z_id, GZ(sz+1:ez-1),istarting,icount)

        istarting(1) = sx+1
        icount(1) = ex-sx-1
        istarting(2) = sy+1
        icount(2) = ey-sy-1
        istarting(3) = sz+1
        icount(3) = ez-sz-1
        r=NF90_put_var(ncdf_id,re_id,DBLE(GRID(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1)),istarting,icount)
        r=NF90_put_var(ncdf_id,im_id,DIMAG(GRID(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1)),istarting,icount)
        r=NF90_put_var(ncdf_id,pot_id,POT(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),istarting,icount)
        call handle_err(r)
    end subroutine

    subroutine write_wf_file_FFTW
        implicit none
        integer :: r,i
        istarting(1) = 1
        icount(1) = NX
        r=NF90_put_var(ncdf_id,x_id,GX,istarting,icount)
        istarting(1) = 1
        icount(1) = NY
        r=NF90_put_var(ncdf_id,y_id,GY,istarting,icount)
        istarting(1) = 1+local_k_offset
        icount(1) = ez
        r=NF90_put_var(ncdf_id,z_id,GZ,istarting,icount)

        istarting(1) = 1
        icount(1) = NX
        istarting(2) = 1
        icount(2) = NY
        istarting(3) = 1+local_k_offset
        icount(3) = ez
        r=NF90_put_var(ncdf_id,re_id,DBLE(GRID),istarting,icount)
        r=NF90_put_var(ncdf_id,im_id,DIMAG(GRID),istarting,icount)
        r=NF90_put_var(ncdf_id,pot_id,POT,istarting,icount)
        call handle_err(r)
    end subroutine

    subroutine close_file()
        implicit none
        integer :: r
        r=NF90_close(ncdf_id)
    end subroutine

    !TODO: FFTW form 
    subroutine read_wf_file(fname)
        implicit none
        integer :: r,rwf_ncid,rwf_re_id,rwf_im_id,rwf_pot_id,file_status,file_age,i
        double precision, dimension(sx:ex,sy:ey,sz:ez) :: realgrid,imaggrid,potgrid
        integer, dimension(13) :: file_stats
        character(len=2048) fname,tmp_fname,latest_fname
        !This next section is a bit weird.
        !It's done this way because I can't find a portable way to ls the working dir in fortran.
        if(RANK .eq. 0) then
            write(6,*) "Reading saved state..."
            if (fname .eq. 'latest') then
                write(6,*) "Finding latest file..."
                file_age = 0
                do i = 0, 999999 !largest filename.
                        write(tmp_fname, '(a,i0.6,a)') 'imag.', i ,".nc"
                        call STAT(tmp_fname, file_stats, file_status)
                        if (file_status == 0) then !if file exists
                            if(file_stats(10) >= file_age) then !if last modification newer than saved
                                file_age = file_stats(10)
                                latest_fname = tmp_fname
                            end if
                        end if
                        write(tmp_fname, '(a,i0.6,a)') 'psi.', i ,".nc"
                        call STAT(tmp_fname, file_stats, file_status)
                        if (file_status == 0) then !if file exists
                            if(file_stats(10) >= file_age) then !if last modification newer than saved
                                file_age = file_stats(10)
                                latest_fname = tmp_fname
                            end if
                        end if
                end do
                if(file_age .eq. 0) then
                    write(6,*) "ERROR: Can't find any previous state! Quitting..."
                    call EXIT(1)
                else 
                    write(6,*) "Using file: ", TRIM(latest_fname), " with modification time: ", CTIME(file_age)
                    fname = latest_fname
                end if
            endif
        end if

        call MPI_Bcast(fname, 2048, MPI_CHARACTER, 0, MPI_COMM, IERR)

        if(RANK .eq. 0) then
            write(6,*) "Opening file: ", TRIM(fname)
        end if

        r = NF90_open_par(fname,IOR(nf90_netcdf4,IOR(NF90_NOWRITE,nf90_MPIIO)),MPI_COMM,MPI_INFO_NULL,rwf_ncid)
        call handle_err(r)

        r = NF90_inq_varid(rwf_ncid, "real",  rwf_re_id)
        r = NF90_inq_varid(rwf_ncid, "imag",  rwf_im_id)
        r = NF90_inq_varid(rwf_ncid,  "pot", rwf_pot_id)

        istarting(1) = sx+1
        icount(1) = ex-sx-1
        istarting(2) = sy+1
        icount(2) = ey-sy-1
        istarting(3) = sz+1
        icount(3) = ez-sz-1

        r = NF90_get_var(rwf_ncid, rwf_re_id, realgrid(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),start=istarting,count=icount)
        r = NF90_get_var(rwf_ncid, rwf_im_id, imaggrid(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),start=istarting,count=icount)
        r = NF90_get_var(rwf_ncid, rwf_pot_id, potgrid(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1),start=istarting,count=icount)

        GRID = realgrid + EYE*imaggrid
        POT = potgrid
        r = NF90_close(rwf_ncid)
        if(RANK .eq. 0) then
            write(6,*) "Reload complete!"
        end if
    end subroutine

    subroutine handle_err(status)
        integer, intent ( in) :: status
     
        if(status /= nf90_noerr) then
          print *, trim(nf90_strerror(status))
          stop "Stopped"
        end if
      end subroutine handle_err

end module