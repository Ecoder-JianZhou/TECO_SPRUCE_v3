program test_select
    implicit none
    integer i, j
    do i = 1,365
        select case(i)
            case (1)
                j = 31
                write(*,*) i, j
            case (31)
                j = 60
                write(*,*) i, j
            case (61)
                j = 91
                write(*,*) i, j
            case (91)
                j=120
                write(*,*) i, j
        end select
        write(*,*) i, j
    end do
end program test_select