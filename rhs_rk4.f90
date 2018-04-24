module rhs_RK4
  use workspace
	contains
    subroutine run_checks_RK4
      implicit none
      if(initialCondType == 2) then
        if(RANK .eq. 0) then
          write(6,'(a)') "Error: Highly non-equilibrium initialCondType is not supported with RK4"
          call finalize_parallel
          call exit(1)
        end if
      end if
    end subroutine
    subroutine RK4_step(rt)
        implicit none
        integer :: rt,i,j,k
        !OpenMP parallelised RK4
        call halo_swap(WS(1)%GRID)
        call RK4_gperhs(WS(1)%GRID, GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = sz+1,ez-1
                do j = sy+1,ey-1
                    do i = sx+1,ex-1
                        GRID_T2(i,j,k) = WS(1)%GRID(i,j,k) + 0.5d0*DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = WS(1)%GRID(i,j,k) + DT*GRID_T1(i,j,k)/6.0d0
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
                        GRID_T2(i,j,k) = WS(1)%GRID(i,j,k) + 0.5d0*DT*GRID_T1(i,j,k)
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
                        GRID_T2(i,j,k) = WS(1)%GRID(i,j,k) + DT*GRID_T1(i,j,k)
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
                        WS(1)%GRID(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/6.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        !Renormalise WF after imaginary time decay
        if(rt .eq. 0) then
          call RK4_renormalise
        end if
    end subroutine

    subroutine RK4_gperhs(gt,kk,rt)
        implicit none
        integer :: rt,i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt, kk
        !Homg GPE
        if (RHSType .eq. 0) then
         !$OMP parallel do private (i,j,k) collapse(3)
           do k = sz+1,ez-1
               do j = sy+1,ey-1
                   do i = sx+1,ex-1
                       kk(i,j,k) = -0.5d0*laplacian(gt,i,j,k)&
                                   + gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&
                                   + POT(i,j,k)*gt(i,j,k) - gt(i,j,k)&
                                   + VELX*EYE*ddx(gt,i,j,k)&
                                   + VELY*EYE*ddy(gt,i,j,k)&
                                   + VELZ*EYE*ddz(gt,i,j,k)&
                                   + OMEGA*EYE*(GX(i)*ddy(gt,i,j,k)-GY(j)*ddx(gt,i,j,k))
                   end do
               end do
           end do
        !$OMP end parallel do
        end if
        if (RHSType .eq. 1) then
        !Trapped GPE
        !$OMP parallel do private (i,j,k) collapse(3)
           do k = sz+1,ez-1
               do j = sy+1,ey-1
                   do i = sx+1,ex-1
                       kk(i,j,k) = -0.5d0*laplacian(gt,i,j,k)&
                                   + harm_osc_C*gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&
                                   + POT(i,j,k)*gt(i,j,k) - harm_osc_mu*gt(i,j,k)&
                                   + OMEGA*EYE*(GX(i)*ddy(gt,i,j,k)-GY(j)*ddx(gt,i,j,k))
                   end do
               end do
           end do
        !$OMP end parallel do
        end if

        if (RHSType .eq. 3) then
        !Quasi-periodic GPE
        GRID_E = CONJG(QP_EX)*gt
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz+1,ez-1
          do j = sy+1,ey-1
            do i = sx+1,ex-1
              kk(i,j,k) = gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))-gt(i,j,k)-QP_EX(i,j,k)*d2dx2(GRID_E,i,j,k)
            end do
          end do
        end do
        !$OMP end parallel do
        GRID_E = CONJG(QP_EY)*gt
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz+1,ez-1
          do j = sy+1,ey-1
            do i = sx+1,ex-1
              kk(i,j,k) = kk(i,j,k)-QP_EY(i,j,k)*d2dy2(GRID_E,i,j,k)
            end do
          end do
        end do
        !$OMP end parallel do
        GRID_E = CONJG(QP_EZ)*gt
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz+1,ez-1
          do j = sy+1,ey-1
            do i = sx+1,ex-1
              kk(i,j,k) = kk(i,j,k)-QP_EZ(i,j,k)*d2dz2(GRID_E,i,j,k)
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
        call MPI_Cart_shift(MPI_COMM_GRID,0,1,rs,rd,IERR)
        call MPI_sendrecv(kk(ex-1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rd,0,&
                          kk(sx,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rs, 0, MPI_COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if

        call MPI_sendrecv(kk(sx+1,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rs,1,&
                          kk(ex,sy:ey,sz:ez),(ey-sy+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rd, 1, MPI_COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if

        !Y dim - tag 2 and 3
        call MPI_Cart_shift(MPI_COMM_GRID,1,1,rs,rd,IERR)
        call MPI_sendrecv(kk(sx:ex,ey-1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rd,2,&
                          kk(sx:ex,sy,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rs, 2, MPI_COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ", IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
        call MPI_sendrecv(kk(sx:ex,sy+1,sz:ez),(ex-sx+1)*(ez-sz+1), MPI_DOUBLE_COMPLEX, rs,3,&
                          kk(sx:ex,ey,sz:ez),(ex-sx+1)*(ez-sz+1),MPI_DOUBLE_COMPLEX, rd, 3, MPI_COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
        !Z dim - tag 4 and 5
        call MPI_Cart_shift(MPI_COMM_GRID,2,1,rs,rd,IERR)
        call MPI_sendrecv(kk(sx:ex,sy:ey,ez-1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE_COMPLEX, rd,4,&
                          kk(sx:ex,sy:ey,sz),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE_COMPLEX, rs, 4, MPI_COMM_GRID,MPISTAT,IERR)
        if(IERR .ne. MPI_SUCCESS) then
          write(6,*) "Error running MPI_sendrecv on Rank: ", RANK, ". Error code: ",IERR
          write(6,*) "Something has gone wrong... Quitting."
          CALL EXIT(1)
        end if
        call MPI_sendrecv(kk(sx:ex,sy:ey,sz+1),(ex-sx+1)*(ey-sy+1), MPI_DOUBLE_COMPLEX, rs,5,&
                          kk(sx:ex,sy:ey,ez),(ex-sx+1)*(ey-sy+1),MPI_DOUBLE_COMPLEX, rd, 5, MPI_COMM_GRID,MPISTAT,IERR)
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
        ddx = (BC(gt,i+1,j,k)-BC(gt,i-1,j,k))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function ddy(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        ddy = (BC(gt,i,j+1,k)-BC(gt,i,j-1,k))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function ddz(gt,i,j,k)
        use params
        implicit none
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        integer :: i,j,k
        ddz = (BC(gt,i,j,k+1)-BC(gt,i,j,k-1))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function d2dx2(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        d2dx2 = (-2.0d0*gt(i,j,k) + BC(gt,i+1,j,k)+BC(gt,i-1,j,k))/(DSPACE**2.0d0)
    end function

    COMPLEX*16 function d2dy2(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        d2dy2 = (-2.0d0*gt(i,j,k) + BC(gt,i,j+1,k)+BC(gt,i,j-1,k))/(DSPACE**2.0d0)
    end function

    COMPLEX*16 function d2dz2(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        d2dz2 = (-2.0d0*gt(i,j,k) + BC(gt,i,j,k+1)+BC(gt,i,j,k-1))/(DSPACE**2.0d0)
    end function

    COMPLEX*16 function laplacian(gt,i,j,k)
        use params
        implicit none
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        integer :: i,j,k
        laplacian = (-6.0d0*gt(i,j,k) + BC(gt,i+1,j,k)+BC(gt,i-1,j,k)&
                                      + BC(gt,i,j+1,k)+BC(gt,i,j-1,k)&
                                      + BC(gt,i,j,k+1)+BC(gt,i,j,k-1))/(DSPACE**2.0d0)
    end function

    complex*16 function BC(gt,i,j,k)
        use params
        implicit none
        complex*16, dimension(sx:ex,sy:ey,sz:ez) :: gt
        integer :: i,j,k,ii,jj,kk
        !Note - Ghost points means that periodic is the default
        ii = i
        jj = j
        kk = k

        !Reflective
        if(i == NX+1 .and. BCX == 0) ii=NX
        if(i == 0   .and. BCX == 0) ii=1
        if(j == NY+1 .and. BCY == 0) jj=NY
        if(j == 0   .and. BCY == 0) jj=1
        if(k == NZ+1 .and. BCZ == 0) kk=NZ
        if(k == 0   .and. BCZ == 0) kk=1

        BC = gt(ii,jj,kk)
        
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

        !Quasi-periodic
        if(i == NX+1 .and. BCX==3) then
          BC = BC*exp(EYE*PI*NVORTZ*GY(jj)/((NY-1)*DSPACE)-EYE*PI*NVORTY*GZ(kk)/((NZ-1)*DSPACE))
        end if
        if(i == 0 .and. BCX==3) then
          BC =BC/exp(EYE*PI*NVORTZ*GY(jj)/((NY-1)*DSPACE)-EYE*PI*NVORTY*GZ(kk)/((NZ-1)*DSPACE))
        end if
        if(j == NY+1 .and. BCY==3) then
          BC = BC*exp(-EYE*PI*NVORTZ*GX(ii)/((NX-1)*DSPACE)+EYE*PI*NVORTX*GZ(kk)/((NZ-1)*DSPACE))
        end if
        if(j == 0 .and. BCY==3) then
          BC = BC/exp(-EYE*PI*NVORTZ*GX(ii)/((NX-1)*DSPACE)+EYE*PI*NVORTX*GZ(kk)/((NZ-1)*DSPACE))
        end if
        if(k == NZ+1 .and. BCZ==3) then
          BC = BC*exp(EYE*PI*NVORTY*GX(ii)/((NX-1)*DSPACE)-EYE*PI*NVORTX*GY(jj)/((NY-1)*DSPACE)&
                                    +EYE*PI*NVORTX*NVORTY)
        end if
        if(k == 0 .and. BCZ==3) then
          BC = BC/exp(EYE*PI*NVORTY*GX(ii)/((NX-1)*DSPACE)-EYE*PI*NVORTX*GY(jj)/((NY-1)*DSPACE)&
                                    +EYE*PI*NVORTX*NVORTY)
        end if
    end function

    subroutine RK4_renormalise
      implicit none
      include 'mpif.h'
      double precision :: local_norm,total_norm,local_pot_int,total_pot_int

      GRID_T1(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1) = WS(1)%GRID(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1)*&
                                                  conjg(WS(1)%GRID(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))

      local_norm=sum(GRID_T1(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))*DSPACE*DSPACE*DSPACE
      call MPI_Allreduce(local_norm, total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW,IERR)

      GRID_T1(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1) = 0.0d0
      where (POT<1.0d0) GRID_T1(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1)=1.0d0

      local_pot_int=sum(GRID_T1(sx+1:ex-1,sy+1:ey-1,sz+1:ez-1))*DSPACE*DSPACE*DSPACE
      call MPI_Allreduce(local_pot_int, total_pot_int, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_FFTW,IERR)

      WS(1)%GRID(sx:ex,sy:ey,sz:ez) = WS(1)%GRID(sx:ex,sy:ey,sz:ez)/sqrt(total_norm)*sqrt(total_pot_int)

    end subroutine
end module