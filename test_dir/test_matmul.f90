program matmultest
    implicit none
    real a(3,3), b(3,1), c(3,1), d(3)
    integer i, j
    do i = 1,3
        b(i,1) = i
        do j = 1,3
            a(i,j) = i+j
        end do
    end do
    d = (/1,2,4/)
    c = matmul(a,b)
    write(*,*) a
    write(*,*) b
    write(*,*) c
    

end program matmultest