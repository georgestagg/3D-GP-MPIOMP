subroutine calc_POT
    use workspace
    implicit none
    integer :: i, j ,k, n
    double precision :: xp,yp,zp

    if (enableTrap) then
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = sz,ez
            do j = sy,ey
                do i = sx,ex
                    POT(i,j,k) = 0.5d0*(TXSCALE*(GX(i)-TX)**2.0d0+TYSCALE*(GY(j)-TY)**2.0d0+TZSCALE*(GZ(k)-TZ)**2.0d0)
                end do
            end do
        end do
        !$OMP end parallel do
    end if
    
end subroutine