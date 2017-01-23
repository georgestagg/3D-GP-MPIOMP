subroutine calc_POT
    use params
    implicit none
    integer :: i, j ,k

    if (enableTrap) then
       !$OMP parallel do private (i,j,k) collapse(3)
        do k = PSZ,PEZ
            do j = PSY,PEY
                do i = PSX,PEX
                    POT(i,j,k) = 0.5d0*(TXSCALE*(GX(i)-TX)**2.0d0+TYSCALE*(GY(j)-TY)**2.0d0+TZSCALE*(GZ(k)-TZ)**2.0d0)
                end do
            end do
        end do
        !$OMP end parallel do
    end if
    
end subroutine