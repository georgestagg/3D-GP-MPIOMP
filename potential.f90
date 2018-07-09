subroutine calc_POT
    use workspace
    implicit none

    POT=0.0d0

    if (enablePot) then
        if(potType .eq. 0) then
            call calc_OBJPOT_obj
        end if
        if(potType .eq. 1) then
            call calc_OBJPOT_surf
        end if
    end if
    if (enableTrap) then
        if (trapType .eq. 0) then
            call add_harmonic_trap
        end if
        if (trapType .eq. 1) then
            call add_hard_circle_trap
        end if
        if (trapType .eq. 2) then
            call add_hard_box_trap
        end if
        if (trapType .eq. 3) then
            call add_soft_circle_trap
        end if
    end if
    
end subroutine

subroutine calc_OBJPOT_surf
    use workspace
    use io
    implicit none
    call read_surf_file_RK4(SURFfilename)
end subroutine

subroutine calc_OBJPOT_obj
    use workspace
    implicit none
    integer :: i,j,k
    double precision :: rx,ry,rz
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                rx = GX(i)-OBJX
                ry = GY(j)-OBJY
                rz = GZ(k)-OBJZ
                POT(i,j,k) = OBJHEIGHT*EXP(-(1.0d0/OBJXSCALE**2.0d0)*(rx**2.0d0) &
                                           -(1.0d0/OBJYSCALE**2.0d0)*(ry**2.0d0) &
                                           -(1.0d0/OBJZSCALE**2.0d0)*(rz**2.0d0))
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine


subroutine add_harmonic_trap
    use workspace
    implicit none
    integer :: i, j ,k
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                POT(i,j,k) = POT(i,j,k) + 0.5d0*(TXSCALE*(GX(i)-TX)**2.0d0+TYSCALE*(GY(j)-TY)**2.0d0+TZSCALE*(GZ(k)-TZ)**2.0d0)
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine

subroutine add_hard_circle_trap
    use workspace
    implicit none
    integer :: i,j,k
    double precision :: rx,ry,rz,r
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                rx = 1.0d0/TXSCALE*(GX(i)-TX)
                ry = 1.0d0/TYSCALE*(GY(j)-TY)
                rz = 1.0d0/TZSCALE*(GZ(k)-TZ)
                r = SQRT(rx**2.0+ry**2.0+rz**2.0)
                if (r > TRAPR) then
                    POT(i,j,k) = POT(i,j,k)+TRAPHEIGHT
                end if
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine

subroutine add_soft_circle_trap
    use workspace
    implicit none
    integer :: i,j,k
    double precision :: rx,ry,rz,r
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                rx = 1.0d0/TXSCALE*(GX(i)-TX)
                ry = 1.0d0/TYSCALE*(GY(j)-TY)
                rz = 1.0d0/TZSCALE*(GZ(k)-TZ)
                r = SQRT(rx**2.0+ry**2.0+rz**2.0)
                POT(i,j,k) = POT(i,j,k)+TRAPHEIGHT/(1.0d0+EXP(-TRAPBETA*(ABS(r)-TRAPR)))
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine

subroutine add_hard_box_trap
    use workspace
    implicit none
    integer :: i,j,k
    double precision :: rx,ry,rz,r
    !$OMP parallel do private (i,j,k) collapse(3)
    do k = sz,ez
        do j = sy,ey
            do i = sx,ex
                rx = 1.0d0/TXSCALE*(GX(i)-TX)
                ry = 1.0d0/TYSCALE*(GY(j)-TY)
                rz = 1.0d0/TZSCALE*(GZ(k)-TZ)
                r = SQRT(rx**2.0+ry**2.0+rz**2.0)
                if (abs(rx) > TRAPR .OR. abs(ry) > TRAPR .OR. abs(rz) > TRAPR) then
                    POT(i,j,k) = POT(i,j,k)+TRAPHEIGHT
                end if
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine