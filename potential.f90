subroutine calc_POT
    use params
    use parallel
    implicit none
    integer :: i, j ,k, n
    double precision :: xp,yp,zp

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

    !many pinning sites
    call srand(rank)
    do n = 1,2000
        xp = RAND()
        yp = RAND()
        zp = RAND()
        
        !$OMP parallel do private (i,j,k) collapse(3)
        do k = PSZ,PEZ
            do j = PSY,PEY
                do i = PSX,PEX
                    POT(i,j,k) = POT(i,j,k) + 100.0d0*EXP((-(GX(i) - (PSX+xp*(PEX-PSX))*DSPACE)**2.0d0 &
                                                           -(GY(j) - (PSY+yp*(PEY-PSY))*DSPACE)**2.0d0 &
                                                           -(GZ(k) - (PSZ+zp*(PEZ-PSZ))*DSPACE)**2.0d0)/0.2d0**2.0d0)
                end do
            end do
        end do
        !$OMP end parallel do
    end do
    
end subroutine