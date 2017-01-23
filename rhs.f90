module rhs
    use params
    use parallel
	contains
    subroutine RK4_step(rt)
        implicit none
        integer :: rt,BC,i,j,k
        !OpenMP parallelised RK4
        call halo_swap(GRID)
        call gperhs(GRID, GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ+1,PEZ-1
                do j = PSY+1,PEY-1
                    do i = PSX+1,PEX-1
                        GRID_T2(i,j,k) = GRID(i,j,k) + 0.5d0*DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = GRID(i,j,k) + DT*GRID_T1(i,j,k)/6.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        call halo_swap(GRID_T2)
        call gperhs(GRID_T2,GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ+1,PEZ-1
                do j = PSY+1,PEY-1
                    do i = PSX+1,PEX-1
                        GRID_T2(i,j,k) = GRID(i,j,k) + 0.5d0*DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/3.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        call halo_swap(GRID_T2)
        call gperhs(GRID_T2,GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ+1,PEZ-1
                do j = PSY+1,PEY-1
                    do i = PSX+1,PEX-1
                        GRID_T2(i,j,k) = GRID(i,j,k) + DT*GRID_T1(i,j,k)
                        GRID_T3(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/3.0d0
                    end do
                end do
            end do
        !$OMP end parallel do

        call halo_swap(GRID_T2)
        call gperhs(GRID_T2,GRID_T1,rt)
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ+1,PEZ-1
                do j = PSY+1,PEY-1
                    do i = PSX+1,PEX-1
                        GRID(i,j,k) = GRID_T3(i,j,k) + DT*GRID_T1(i,j,k)/6.0d0
                    end do
                end do
            end do
        !$OMP end parallel do
    end subroutine

    subroutine gperhs(gt,kk,rt)
        implicit none
        integer :: rt,i,j,k
        complex*16, dimension(PSX:PEX,PSY:PEY,PSZ:PEZ) :: gt, kk

        !OpenMP parallelised Homg GPE
       !!$OMP parallel do private (i,j,k) collapse(3)
       !    do k = PSZ+1,PEZ-1
       !        do j = PSY+1,PEY-1
       !            do i = PSX+1,PEX-1
       !                kk(i,j,k) = -0.5d0*(-6.0d0*gt(i,j,k)&
       !                            + gt(BC(i+1,0),j,k)+gt(BC(i-1,0),j,k)&
       !                            + gt(i,BC(j+1,1),k)+gt(i,BC(j-1,1),k)&
       !                            + gt(i,j,BC(k+1,2))+gt(i,j,BC(k-1,2)))/(DSPACE**2.0d0)&
       !                            + gt(i,j,k)*gt(i,j,k)*CONJG(gt(i,j,k))&
       !                            + POT(i,j,k)*gt(i,j,k) - gt(i,j,k)&
       !                            + OMEGA*EYE*((GX(i)-(NX-1)*DSPACE/2.0d0)*ddy(gt,i,j,k)&
       !                            - (GY(j)-(NY-1)*DSPACE/2.0d0)*ddx(gt,i,j,k))
       !            end do
       !        end do
       !    end do
       !!$OMP end parallel do

       !$OMP parallel do private (i,j,k) collapse(3)
           do k = PSZ+1,PEZ-1
               do j = PSY+1,PEY-1
                   do i = PSX+1,PEX-1
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

        !Damping
        if(dble(GAMMAC) > 0.0d0) then
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ+1,PEZ-1
                do j = PSY+1,PEY-1
                    do i = PSX+1,PEX-1
                        kk(i,j,k)=kk(i,j,k)/(EYE-GAMMAC)
                    end do
                end do
            end do
        !$OMP end parallel do
        else
        !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ+1,PEZ-1
                do j = PSY+1,PEY-1
                    do i = PSX+1,PEX-1
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
        complex*16, dimension(PSX:PEX,PSY:PEY,PSZ:PEZ) :: kk
        include 'mpif.h'

        !X dim - tag 0 and 1
        call MPI_Cart_shift(COMM_GRID,0,1,rs,rd,IERR)
        call MPI_sendrecv(kk(PEX-1,PSY:PEY,PSZ:PEZ),(PEY-PSY+1)*(PEZ-PSZ+1), MPI_DOUBLE_COMPLEX, rd,0,&
                          kk(PSX,PSY:PEY,PSZ:PEZ),(PEY-PSY+1)*(PEZ-PSZ+1),MPI_DOUBLE_COMPLEX, rs, 0, COMM_GRID,MPISTAT,IERR)
        call MPI_sendrecv(kk(PSX+1,PSY:PEY,PSZ:PEZ),(PEY-PSY+1)*(PEZ-PSZ+1), MPI_DOUBLE_COMPLEX, rs,1,&
                          kk(PEX,PSY:PEY,PSZ:PEZ),(PEY-PSY+1)*(PEZ-PSZ+1),MPI_DOUBLE_COMPLEX, rd, 1, COMM_GRID,MPISTAT,IERR)
        !Y dim - tag 2 and 3
        call MPI_Cart_shift(COMM_GRID,1,1,rs,rd,IERR)
        call MPI_sendrecv(kk(PSX:PEX,PEY-1,PSZ:PEZ),(PEX-PSX+1)*(PEZ-PSZ+1), MPI_DOUBLE_COMPLEX, rd,2,&
                          kk(PSX:PEX,PSY,PSZ:PEZ),(PEX-PSX+1)*(PEZ-PSZ+1),MPI_DOUBLE_COMPLEX, rs, 2, COMM_GRID,MPISTAT,IERR)
        call MPI_sendrecv(kk(PSX:PEX,PSY+1,PSZ:PEZ),(PEX-PSX+1)*(PEZ-PSZ+1), MPI_DOUBLE_COMPLEX, rs,3,&
                          kk(PSX:PEX,PEY,PSZ:PEZ),(PEX-PSX+1)*(PEZ-PSZ+1),MPI_DOUBLE_COMPLEX, rd, 3, COMM_GRID,MPISTAT,IERR)
        !Z dim - tag 4 and 5
        call MPI_Cart_shift(COMM_GRID,2,1,rs,rd,IERR)
        call MPI_sendrecv(kk(PSX:PEX,PSY:PEY,PEZ-1),(PEX-PSX+1)*(PEY-PSY+1), MPI_DOUBLE_COMPLEX, rd,4,&
                          kk(PSX:PEX,PSY:PEY,PSZ),(PEX-PSX+1)*(PEY-PSY+1),MPI_DOUBLE_COMPLEX, rs, 4, COMM_GRID,MPISTAT,IERR)
        call MPI_sendrecv(kk(PSX:PEX,PSY:PEY,PSZ+1),(PEX-PSX+1)*(PEY-PSY+1), MPI_DOUBLE_COMPLEX, rs,5,&
                          kk(PSX:PEX,PSY:PEY,PEZ),(PEX-PSX+1)*(PEY-PSY+1),MPI_DOUBLE_COMPLEX, rd, 5, COMM_GRID,MPISTAT,IERR)
    end subroutine

    COMPLEX*16 function ddx(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(PSX:PEX,PSY:PEY,PSZ:PEZ) :: gt
        ddx = (gt(BC(i+1,0),j,k)-gt(BC(i-1,0),j,k))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function ddy(gt,i,j,k)
        use params
        implicit none
        integer :: i,j,k
        complex*16, dimension(PSX:PEX,PSY:PEY,PSZ:PEZ) :: gt
        ddy = (gt(i,BC(j+1,0),k)-gt(i,BC(j-1,0),k))/(2.0d0*DSPACE)
    end function

    COMPLEX*16 function ddz(gt,i,j,k)
        use params
        implicit none
        complex*16, dimension(PSX:PEX,PSY:PEY,PSZ:PEZ) :: gt
        integer :: i,j,k
        ddz = (gt(i,j,BC(k+1,0))-gt(i,j,BC(k-1,0)))/(2.0d0*DSPACE)
    end function

    integer function BC(s,n)
        use params
        implicit none
        integer :: s,n
        BC=s
        select case (n)
            case (0)
                if(s == PEX+1) BC=PSX
                if(s == PSX-1) BC=PEX
            case (1)
                if(s == PEY+1) BC=PSY
                if(s == PSY-1) BC=PEY
            case (2)
                if(s == PEZ+1) BC=PSZ
                if(s == PSZ-1) BC=PEZ
        end select
    end function
end module