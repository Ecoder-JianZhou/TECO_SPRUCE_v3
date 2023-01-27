program cycle
    implicit none
    integer j
    
    do j=1,5
        call test()
    enddo

end program cycle

subroutine test()
    implicit none
    real, save :: a = 0.
    integer i

    do i=1,5
        a = a+1
    enddo
    write(*,*) "a = ", a
end subroutine test