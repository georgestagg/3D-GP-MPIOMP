module rhs_RK4
    use workspace
	contains
    subroutine RK4_step(rt)
        implicit none
        integer :: rt,BC,i,j,k
        !OpenMP parallelised RK4
        call halo_swap(GRID)
        call RK4_gperhs(GRID, GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        GRID_T2(i,j,k) = GRID(i,j,k) + 0.5d0*DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = GRID(i,j,k) + DT*GRID_T1(i,j,k)/6.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        call halo_swap(GRID_T2)
        call RK4_gperhs(GRID_T2,GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        GRID_T2(i,j,k) = GRID(i,j,k) + 0.5d0*DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/3.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        call halo_swap(GRID_T2)
        call RK4_gperhs(GRID_T2,GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        GRID_T2(i,j,k) = GRID(i,j,k) + DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/3.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        call halo_swap(GRID_T2)
        call RK4_gperhs(GRID_T2,GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        GRID(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/6.0d0
                    end do
                end do
            end do
        !$OMP end parallel do
    end subroutine

    subroutine RK4_gperhs(gt,kk,rt)
        implicit none
        integer :: rt,i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt, kk
        !OpenMP parallelised Homg GPE
        if (RHSType .eq. 0) then
         !$OMP parallel do private (i,j,k) collapse(3)
           do k = sz+1,ez-1
               do j = sy+1,ey-1
                   do i = sx+1,ex-1
                       kk(i,j,k) = -0.5d0*(-6.0d0*gt(i,j,k)&
                                   + gt(BC(i+1,0),j,k)+gt(BC(i-1,0),j,k)&
                                   + gt(i,BC(j+1,1),k)+gt(i,BC(j-1,1),k)&
                                   + gt(i,j,BC(k+1,2))+gt(i,j,BC(k-1,2)))/(DSPACE**2.0d0)&
                                   + gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&
                                   + POT(i,j,k)*gt(i,j,k) - gt(i,j,k)&
                                   + VELX*EYE*ddx(gt,i,j,k)&
                                   + VELY*EYE*ddy(gt,i,j,k)&
                                   + VELZ*EYE*ddz(gt,i,j,k)&
                                   + OMEGA*EYE*((GX(i)-(NX-1)*DSPACE/2.0d0)*ddy(gt,i,j,k)&
                                   - (GY(j)-(NY-1)*DSPACE/2.0d0)*ddx(gt,i,j,k))
                   end do
               end do
           end do
        !$OMP end parallel do
        end if
        if (RHSType .eq. 1) then
          !OpenMP parallelised trapped GPE
        !$OMP parallel do private (i,j,k) collapse(3)
           do k = sz+1,ez-1
               do j = sy+1,ey-1
                   do i = sx+1,ex-1
                       kk(i,j,k) = -0.5d0*(-6.0d0*gt(i,j,k)&
                                   + gt(BC(i+1,0),j,k)+gt(BC(i-1,0),j,k)&
                                   + gt(i,BC(j+1,1),k)+gt(i,BC(j-1,1),k)&
                                   + gt(i,j,BC(k+1,2))+gt(i,j,BC(k-1,2)))/(DSPACE**2.0d0)&
                                   + harm_osc_C*gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&
                                   + POT(i,j,k)*gt(i,j,k) - harm_osc_mu*gt(i,j,k)&
                                   + OMEGA*EYE*((GX(i)-(NX-1)*DSPACE/2.0d0)*ddy(gt,i,j,k)&
                                   - (GY(j)-(NY-1)*DSPACE/2.0d0)*ddx(gt,i,j,k))
                   end do
               end do
           end do
        !$OMP end parallel do
        end if

        !Damping
        if(dble(GAMMAC) > 0.0d0) then
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        kk(i,j,k)=kk(i,j,k)/(EYE-GAMMAC)
                    end do
                end do
            end do
        !$OMP end parallel do
        else
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        kk(i,j,k) = kk(i,j,k)/EYE
                    end do
                end do
            end do
        !$OMP end parallel do
        end if
    end subroutine

    subroutine halo_swap(kk)
        use parallel
        implicit none
        integer :: rs,rd
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: kk
        include 'mpif.h'

        !X dim - tag 0 and 1
        call MPI_Cart_shift(COMM_GRID,0,1,rs,rd,IERR)
        call MPI_sendrecv(kk(ex-1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rd,0,&
                          kk(sx,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rs, 0, COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if

        call MPI_sendrecv(kk(sx+1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rs,1,&
                          kk(ex,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rd, 1, COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if

        !Y dim - tag 2 and 3
        call MPI_Cart_shift(COMM_GRID,1,1,rs,rd,IERR)
        call MPI_sendrecv(kk(sx:ex,ey-1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rd,2,&
                          kk(sx:ex,sy,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rs, 2, COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
        call MPI_sendrecv(kk(sx:ex,sy+1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rs,3,&
                          kk(sx:ex,ey,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rd, 3, COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
        !Z dim - tag 4 and 5
        call MPI_Cart_shift(COMM_GRID,2,1,rs,rd,IERR)
        call MPI_sendrecv(kk(sx:ex,sy:ey,ez-1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE_COMPLEX, rd,4,&
                          kk(sx:ex,sy:ey,sz),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE_COMPLEX, rs, 4, COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
        call MPI_sendrecv(kk(sx:ex,sy:ey,sz+1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE_COMPLEX, rs,5,&
                          kk(sx:ex,sy:ey,ez),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE_COMPLEX, rd, 5, COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
    end subroutine

    COMPLEX*16 function ddx(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        ddx = (gt(BC(i+1,0),j,k)-gt(BC(i-1,0),j,k))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function ddy(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        ddy = (gt(i,BC(j+1,1),k)-gt(i,BC(j-1,1),k))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function ddz(gt,i,j,k)
        use params
        implicit none
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        integer :: i,j,k
        ddz = (gt(i,j,BC(k+1,2))-gt(i,j,BC(k-1,2)))/(2.0d0*DSPACE)
    end function

    integer function BC(s,n)
        use params
        implicit none
        integer :: s,n
        BC=s
        select case (n)
            case (0)
                if(s == ex+1) BC=sx
                if(s == sx-1) BC=ex
            case (1)
                if(s == ey+1) BC=sy
                if(s == sy-1) BC=ey
            case (2)
                if(s == ez+1) BC=sz
                if(s == sz-1) BC=ez
        end select
    end function
end module