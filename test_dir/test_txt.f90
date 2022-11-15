program readtxt
    implicit none
    integer COUNT, STAT, n
    character(12) a
    character b,c,d,e
    real g(4)
    COUNT=0
    OPEN(1,FILE="filelist.txt")
    DO WHILE (.TRUE.)
        READ(1,*,IOSTAT=STAT)(g(n), n=1,4)
        IF(STAT .NE. 0) EXIT
        COUNT=COUNT+1
        write(*,*) COUNT, g
    ENDDO
    write(*,*) COUNT
    CLOSE(1)
    write(*,*) COUNT
end program readtxt