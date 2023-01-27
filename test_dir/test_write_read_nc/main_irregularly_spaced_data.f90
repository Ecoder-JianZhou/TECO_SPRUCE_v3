program ex
    implicit none
    real(kind=4), DIMENSION(:,:), ALLOCATABLE :: positions
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: mydata
    REAL(KIND=4) :: x, y
    INTEGER(KIND=4) :: NX, NY, naxes
    INTEGER(KIND=4) :: J
    CHARACTER(LEN=50) :: outfile

    outfile = 'ex12.nc'
    NX = 2
    NY = 2
    naxes = 2
    ALLOCATE(mydata(NX*NY))
    ALLOCATE(positions(naxes,NX*NY))
    positions(1:2,1) = (/0.0, 0.0/)
    positions(1:2,2) = (/0.0, 0.2/)
    positions(1:2,3) = (/0.3, 0.5/)
    positions(1:2,4) = (/0.4, 0.9/)
    DO J=1,NX*NY
        x = positions(1,J)
        y = positions(2,J)
        mydata(J) = cos(sqrt(x*x+y*y)/10.0)
    ENDDO
    !write netCDF file
    CALL wirrgrid(outfile,positions,mydata,NX*NY,naxes)

end program ex

!WIRRGRID - write a netCDF gridfile
!:=========================================================================
SUBROUTINE wirrgrid(outfile,positions,idata,pointnums,axes)
    USE netcdf
    IMPLICIT NONE
    REAL(KIND=4), DIMENSION(axes,pointnums), INTENT(IN) :: positions
    REAL(KIND=4), DIMENSION(pointnums), INTENT(IN) :: idata
    INTEGER(KIND=4) :: ncid, x_dimid, y_dimid, d_dimid
    INTEGER(KIND=4) :: x_varid, y_varid, varid
    INTEGER(KIND=4), DIMENSION(2) :: dimids
    INTEGER(KIND=4), INTENT(IN) :: pointnums, axes
    INTEGER(KIND=4) :: J, K, c
    CHARACTER(LEN=50), INTENT(IN) :: outfile

    !Create the netCDF file.
    CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))
    !Define the dimensions.
    CALL check(nf90_def_dim(ncid, "pointnums", pointnums, x_dimid))
    CALL check(nf90_def_dim(ncid, "axes", axes, y_dimid))
    !Dimension ID's
    dimids = (/ y_dimid, x_dimid /)
    !Define coordinate variables
    CALL check(nf90_def_var(ncid, "grid", NF90_FLOAT, dimids, x_varid))
    !Define data variable
    CALL check(nf90_def_var(ncid, "Perturbations", NF90_FLOAT, x_dimid, varid))
    !Add attributes
    CALL check(nf90_put_att(ncid,varid,"units","%"))
    CALL check(nf90_put_att(ncid,varid,"field","Perturbations, vector"))
    CALL check(nf90_put_att(ncid,varid,"positions","grid"))
    CALL check(nf90_enddef(ncid)) !End Definitions
    !Write Data
    CALL check(nf90_put_var(ncid, x_varid, positions))
    CALL check(nf90_put_var(ncid, varid, idata))
    CALL check(nf90_close(ncid))
END SUBROUTINE wirrgrid
!:=========================================================================

! check (ever so slightly modified from www.unidata.ucar.edu)
subroutine check(istatus)
    use netcdf
    implicit none
    integer, intent(in) :: istatus
    if(istatus /= nf90_noerr) then
        write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
end subroutine check