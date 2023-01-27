program test
    use example
    implicit none
    integer(kind=4) :: NX, NY
    character(len=50) :: infile
    real, allocatable :: indata(:,:), xpos(:), ypos(:)
    

    call sub1()
    infile = "test.nc"
    call griddims(infile, NX, NY)
    write(*,*) NX, NY
    allocate(indata(NX,NY))
    allocate(xpos(NX))
    allocate(ypos(NY))
    call readgrid(infile, xpos, ypos, indata, NX, NY)
    write(*,*) indata
    deallocate(indata)
    deallocate(xpos)
    deallocate(ypos)
end program test