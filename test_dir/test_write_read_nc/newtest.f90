module example
    contains
    subroutine sub1()
        use netcdf
        implicit none
        write(*,*) "example-sub1 ..."
    end subroutine sub1

    ! GRIDDIMS - Get dimensions of a netCDF 2D gridfile
    !===================================================
    subroutine griddims(infile, NX, NY)
        use netcdf
        implicit none
        integer(kind=4), intent(out) :: NX, NY
        integer(kind=4) :: ncid
        character(len=50), intent(in) :: infile
        character(len=50) :: xname, yname
        ! Open netCDF file
        call check(nf90_open(infile, nf90_nowrite, ncid))

        ! Inquire about the dimensions
        ! :--------:--------:--------:--------:--------:--------
        call check(nf90_inquire_dimension(ncid,1,xname, NX))
        call check(nf90_inquire_dimension(ncid,2,yname, NY))

        ! close netCDF file
        call check(nf90_close(ncid))
    end subroutine griddims

    ! READGRID - read a netCDF gridfile
    !: ====================================================
    subroutine readgrid(infile, xpos, ypos, idata, NX, NY)
        use netcdf
        implicit none
        real(kind=4), dimension(NX), intent(out) :: xpos
        real(kind=4), dimension(NY), intent(out) :: ypos
        real(kind=4), dimension(NX, NY), intent(out) :: idata
        integer(kind=4), intent(in) :: NX, NY
        integer(kind=4), dimension(2) :: dimids
        integer(kind=4) :: ncid, xtype, ndims, varid
        character(len=50), intent(in) :: infile
        character(len=50) :: xname, yname, vname

        ! Open netCDF file
        ! :--------:--------:--------:--------:--------:--------
        call check(nf90_open(infile, nf90_nowrite, ncid))

        ! Get the values of the coordinates and put them in xpos and ypos
        ! :--------:--------:--------:--------:--------:--------
        call check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
        call check(nf90_inq_varid(ncid, vname, varid))
        call check(nf90_get_var(ncid, varid, xpos))

        call check(nf90_inquire_variable(ncid,2,vname, xtype, ndims, dimids))
        call check(nf90_inq_varid(ncid, vname, varid))
        call check(nf90_get_var(ncid, varid, ypos))

        ! Get the values of the perturbations and put them in idata
        ! :--------:--------:--------:--------:--------:--------
        call check(nf90_inquire_variable(ncid,2,vname,xtype,ndims,dimids))
        call check(nf90_inq_varid(ncid, vname, varid))
        call check(nf90_get_var(ncid, varid, idata))

        ! Close netCDF file
        ! :--------:--------:--------:--------:--------:--------
        call check(nf90_close(ncid))
    end subroutine readgrid

    ! check (ever so slightly modified from www.unidata.ucar.edu)
    subroutine check(istatus)
        use netcdf
        implicit none
        integer, intent(in) :: istatus
        if(istatus /= nf90_noerr) then
            write(*,*) trim(adjustl(nf90_strerror(istatus)))
        end if
    end subroutine check
end module example