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
    if (enablePins) then
        call SRAND(randSeedPins)
        do n = 1,numPotPins
            xp = RAND()
            yp = RAND()
            zp = RAND()
            
            !$OMP parallel do private (i,j,k) collapse(3)
            do k = PSZ,PEZ
                do j = PSY,PEY
                    do i = PSX,PEX
                        POT(i,j,k) = POT(i,j,k) + pinPotDepth*EXP(( -(GX(i) - (xp*NX)*DSPACE)**2.0d0 &
                                                                    -(GY(j) - (yp*NY)*DSPACE)**2.0d0 &
                                                                    -(GZ(k) - (zp*NZ)*DSPACE)**2.0d0)/pinSigma**2.0d0)
                    end do
                end do
            end do
            !$OMP end parallel do
        end do
    end if
    
end subroutine