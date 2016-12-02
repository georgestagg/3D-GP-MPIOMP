module output
    use params
    use parallel
    use netcdf
    implicit none
    include 'mpif.h'
    integer :: ncdf_id,x_id,y_id,z_id,re_id,im_id,pot_id
    integer       :: icount(3)
    integer       :: istarting(3)
    Contains
    subroutine make_file(fname)
        implicit none
        integer :: dims(3)
        integer :: r,x_dim_id,y_dim_id,z_dim_id
        character(len=*) :: fname
        r=NF90_create(fname , IOR(nf90_netcdf4, nf90_classic_model), ncdf_id, comm=COMM_GRID, info=MPI_INFO_NULL)
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

        r=NF90_enddef(ncdf_id)
        r=NF90_sync(ncdf_id)

    end subroutine

    subroutine write_wf_file()
        implicit none
        integer :: r,i
        r=nf90_var_par_access(ncdf_id, x_id, NF90_COLLECTIVE)
        !Dont forget to ignore ghost points!
        istarting(1) = PSX+1
        istarting(2) = PSY+1
        istarting(3) = PSZ+1
        icount(1) = PEX-PSX-2
        r=NF90_put_var(ncdf_id, x_id, GX(PSX+1:PEX-1),istarting,icount)
    end subroutine

    subroutine close_file()
        implicit none
        integer :: r
        r=NF90_close(ncdf_id)
    end subroutine

    subroutine handle_err(status)
        integer, intent ( in) :: status
     
        if(status /= nf90_noerr) then
          print *, trim(nf90_strerror(status))
          stop "Stopped"
        end if
      end subroutine handle_err

end module