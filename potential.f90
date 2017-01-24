subroutine calc_POT
    use params
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


    do n = 1,10000
        xp = RAND()
        yp = RAND()
        zp = RAND()

        do k = PSZ,PEZ
            do j = PSY,PEY
                do i = PSX,PEX
                    POT(i,j,k) = POT(i,j,k) + 100.0d0*EXP((-(GX(i) - xp*(NX-1)*DSPACE)**2.0d0 &
                                                          -(GY(j) - yp*(NY-1)*DSPACE)**2.0d0 &
                                                          -(GZ(k) - zp*(NZ-1)*DSPACE)**2.0d0)/0.1d0**2.0d0)
                end do
            end do
        end do
    end do
    
end subroutine