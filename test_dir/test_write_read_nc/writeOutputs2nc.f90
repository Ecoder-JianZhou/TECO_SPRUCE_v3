module ncdata_io
    implicit none
    
    contains
    ! ----------------------------------------------
    subroutine read_nc()
        use typeSizes
        use netcdf90
        implicit none
        include 'netcdf.inc'
        ! Variables
        character(len=50) filename
        integer ncid,ierr
        ! Body of Console3
        filename ='G:/si_monthly_1850-2012_r2.5x2.5.nc'
        print *, 'Hello World'

        ierr = nf90_open(filename,nf_nowrite,ncid)
        print *, ncid
    end subroutine read_nc

    subroutine write_nc()
      ! This program provides an elementary check of some of the parts of the 
      !   Fortran 90 interface to netCDF 3.5. It is a Fortran 90 implementation
      !   of the nctst.cpp program provided with the C++ interface to netcdf
      !   (in the src/cxx directory of the netcdf distribution). 
        use typeSizes
        use netcdf90
        implicit none
        include 'netcdf.inc'
        
        ! netcdf related variables
        integer :: ncFileID,                                   &
                  latDimID, lonDimID, frTimeDimID, timeDimID, &
                  pressVarID, latVarID, lonVarID, frTimeVarID, refTimeVarID, scalarVarID
                  
        ! Local variables
        integer, parameter :: numLats = 4, numLons = 3, &
                              numFrTimes = 2, timeStringLen = 20
        character (len = *), parameter :: fileName = "example.nc"
        integer :: counter                      
        real, dimension(numLons, numLats, numFrTimes) :: pressure
        integer (kind = FourByteInt), dimension(numFrTimes) :: frTimeVals
        real (kind = FourByteReal) fillVal;
        real (kind = FourByteReal), dimension(2) :: validRange;
        
        ! --------------------
        ! Code begins
        ! --------------------
        if(.not. byteSizesOK()) then
          print *, "Compiler does not appear to support required kinds of variables."
          stop
        end if
          
        ! Create the file
        call check(nf90_create(path = trim(fileName), cmode = nf90_clobber, ncid = ncFileID))
        
        ! Define the dimensions
        call check(nf90_def_dim(ncid = ncFileID, name = "lat",     len = numLats,        dimid = latDimID))
        call check(nf90_def_dim(ncid = ncFileID, name = "lon",     len = numLons,        dimid = lonDimID))
        call check(nf90_def_dim(ncid = ncFileID, name = "frtime",  len = nf90_unlimited, dimid = frTimeDimID))
        call check(nf90_def_dim(ncid = ncFileID, name = "timelen", len = timeStringLen,  dimid = timeDimID))
        ! Create variables and attributes
        call check(nf90_def_var(ncid = ncFileID, name = "P", xtype = nf90_float,     &
                          dimids = (/ lonDimID, latDimID, frTimeDimID /), varID = pressVarID) )
        call check(nf90_put_att(ncFileID, pressVarID, "long_name",   "pressure at maximum wind"))
        call check(nf90_put_att(ncFileID, pressVarID, "units",       "hectopascals") )
        ! Use 4-byte reals explicitly, to match 4-byte attribute type in test file
        validRange(1) = 0.
        validRange(2) = 1500
        call check(nf90_put_att(ncFileID, pressVarID, "valid_range", validRange))
        ! Use a 4-byte float constant, to match variable type
        fillVal = -9999.0
        call check(nf90_put_att(ncFileID, pressVarID,  "_FillValue", fillVal ) )
                            
        call check(nf90_def_var(ncFileID, "lat", nf90_float, dimids = latDimID, varID = latVarID) )
        call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
        call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))
        call check(nf90_def_var(ncFileID, "lon", nf90_float, lonDimID, lonVarID) )
        call check(nf90_put_att(ncFileID, lonVarID, "long_name", "longitude"))
        call check(nf90_put_att(ncFileID, lonVarID, "units",     "degrees_east"))
        call check(nf90_def_var(ncFileID, "frtime", nf90_int, frTimeDimID, frTimeVarID) )
        call check(nf90_put_att(ncFileID, frTimeVarID, "long_name", "forecast time"))
        call check(nf90_put_att(ncFileID, frTimeVarID, "units",     "hours"))
        call check(nf90_def_var(ncFileID, "reftime", nf90_char, timeDimID, refTimeVarID) )
        call check(nf90_put_att(ncFileID, refTimeVarID, "long_name", "reference time"))
        call check(nf90_put_att(ncFileID, refTimeVarID, "units",     "text_time"))
                          
        ! In the C++ interface the define a scalar variable - do we know how to do this? 
        call check(nf90_def_var(ncFileID, "ScalarVariable", nf90_real, scalarVarID))
        
        ! Global attributes
        call check(nf90_put_att(ncFileID, nf90_global, "history", &
                          "created by Unidata LDM from NPS broadcast"))
        call check(nf90_put_att(ncFileID, nf90_global, "title", &
                          "NMC Global Product Set: Pressure at Maximum Wind"))
        
        ! Leave define mode
        call check(nf90_enddef(ncfileID))
        
        ! Write the dimension variables
        call check(nf90_put_var(ncFileID, latVarId,     (/ -90., -87.5, -85., -82.5 /)) )
        call check(nf90_put_var(ncFileID, lonVarId,     (/ -180, -175, -170 /)      ) )
        ! Don't use anonymous array here, in case platform has 8-byte integers
        frTimeVals(1) = 12
        frTimeVals(2) = 18
        call check(nf90_put_var(ncFileID, frTimeVarId,  frTimeVals                  ) )
        call check(nf90_put_var(ncFileID, reftimeVarID, "1992-3-21 12:00"           ) )
        
        ! Write the pressure variable. Write a slab at a time to check incrementing.
        pressure = 949. + real(reshape( (/ (counter, counter = 1, numLats * numLons * numFrTimes) /),  &
                                          (/ numLons, numLats, numFrTimes /) ) )
        call check(nf90_put_var(ncFileID, pressVarID, pressure(:, :, 1:1)) )
        call check(nf90_put_var(ncFileID, pressVarID, pressure(:, :, 2:2), start = (/ 1, 1, 2 /)) )
        
        call check(nf90_put_var(ncFileID, scalarVarID, 10))
        call check(nf90_close(ncFileID))
    end subroutine write_nc


    ! Internal subroutine - checks error status after each netcdf, prints out text message each time
        !   an error code is returned. 
    subroutine check(status)
      integer, intent ( in) :: status
      
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
      end if
    end subroutine check 

end module ncdata_io