program ex
    implicit none
    real(kind=4), dimension(:,:), allocatable :: mydata
    real(kind=4), dimension(:), allocatable :: xpos, ypos
    real(kind=4) :: x, y
    integer(kind=4) :: NX
    integer(kind=4) :: NY
    integer(kind=4) :: J, K
    character(LEN=50) :: outfile

    outfile = "ex1.nc"
    NX = 360
    NY = 180

    allocate(mydata(NX, NY))
    allocate(xpos(NX))
    allocate(ypos(NY))

    ! Populate X and Y grid locations
    x = -179.0
    do J = 1, NX
        xpos(J) = x
        x = x + 1.0
    end do

    y = -90.0
    do J = 1, NY
        ypos(J) = y
        y = y + 1.0
    end do

    ! make up some data
    do J= 1, NX
        x = xpos(J)
        do K = 1, NY
            y = ypos(K)
            mydata(J,K) = cos(sqrt(x*x+y*y)/10)
        end do
    end do

    ! write netcdf file
    call writegrid(outfile, xpos, ypos, mydata, NX, NY)
end program ex

! WRIATEGRID - write a netCDF gridfile
! =========================================================
subroutine writegrid(outfile, xpos, ypos, idata, NX, NY)
    use netcdf
    implicit none
    real(kind=4), dimension(NX), intent(in) :: xpos
    real(kind=4), dimension(NY), intent(in) :: ypos
    real(kind=4), dimension(NX, NY), intent(in) :: idata
    integer(kind=4) :: ncid, x_dimid, y_dimid
    integer(kind=4) :: x_varid, y_varid, varid
    integer(kind=4), dimension(2) :: dimids
    integer(kind=4), intent(in) :: NX, NY
    character(len=50), intent(in) :: outfile

    ! create the netCDF file
    call check(nf90_create(outfile, NF90_CLOBBER, ncid))

    ! define the dimensions.
    call check(nf90_def_dim(ncid, "lon", NX, x_dimid))
    call check(nf90_def_dim(ncid, "lat", NY, y_dimid))

    ! define coordinate variables
    call check(nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, x_varid))
    call check(nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, y_varid))

    ! add units for the coordinate variables
    call check(nf90_put_att(ncid, x_varid, "units", "deg"))
    call check(nf90_put_att(ncid, y_varid, "units", "deg"))

    ! define our data array
    dimids = (/x_dimid, y_dimid/)

    ! define variable
    call check(nf90_def_var(ncid, "Perturbations", NF90_FLOAT, dimids, varid))
    call check(nf90_put_att(ncid, varid, "units", "dVs %"))
    call check(nf90_put_att(ncid, varid, "Function_Type", "Matern"))

    call check(nf90_put_att(ncid, nf90_global, "title", "Realization by KL Expansion"))

    ! end definitions
    call check(nf90_enddef(ncid))  

    ! write data
    call check(nf90_put_var(ncid, x_varid, xpos))
    call check(nf90_put_var(ncid, y_varid, ypos))
    call check(nf90_put_var(ncid, varid, idata))
    call check(nf90_close(ncid))
end subroutine writegrid

! check (ever so slightly modified from www.unidata.ucar.edu)
subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent(in) :: istatus
    if(istatus /= nf90_noerr) then
        write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
end subroutine check