subroutine setupGXYZ
    use params
    implicit none
    integer :: i, j ,k
    do i = PSX,PEX
        GX(i) = (i-1)*DSPACE
    end do
    
    do j = PSY,PEY
        GY(j) = (j-1)*DSPACE
    end do
    
    do k = PSZ,PEZ
        GZ(k) = (k-1)*DSPACE
    end do
end subroutine

subroutine setupMeshgrid
    use params
    implicit none
    integer :: i, j ,k
    do k = PSZ,PEZ
        do j = PSY,PEY
            do i = PSX,PEX
                MGX(i,j,k) = GX(i)
                MGY(i,j,k) = GY(j)
                MGZ(i,j,k) = GZ(k)
            end do
        end do
    end do
end subroutine

subroutine initCond
    use params
    implicit none
    GRID = 1.00d0
end subroutine

subroutine logger(text)
    character(len=*),intent(in)   :: text
    character(len=80) :: fname
    write(fname, '(i0.4,a)') 'output.log'
    open (8, FILE = fname, ACCESS = 'APPEND')
        write (unit=8,fmt="(a)") text
    close(8)
end subroutine