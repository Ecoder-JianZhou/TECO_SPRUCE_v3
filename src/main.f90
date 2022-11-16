program TECO
    use mod_data
    use mod_spinup
    use mod_mcmc
    use driver
    ! to run TECO simulation, spin-up and data simulation
    implicit none

    parafile      = "input/parameters.txt"
    climatefile   = "input/forcing.txt"
    snowdepthfile = "input/SPRUCE_Snow_Depth_2011-2014.txt"
    call get_params()                           ! read parameters values
    call get_forcingdata()                      ! read forcing data
    call get_snowdepth()
    call initialize()                           ! initializations
    if (do_spinup .eq. 1) call run_spinup() 
    if (do_mcmc .eq. 1) then
        call read_obs()
        call run_mcmc()
    endif
    write(*,*)Dheat,gddonset
    call teco_simu()
    ! end of the simulation, then deallocate the forcing_data
    deallocate(forcing%year)
    deallocate(forcing%doy)
    deallocate(forcing%hour)
    deallocate(forcing%Tair)
    deallocate(forcing%Tsoil)
    deallocate(forcing%RH)
    deallocate(forcing%VPD)
    deallocate(forcing%Rain)
    deallocate(forcing%WS)
    deallocate(forcing%PAR)
    deallocate(forcing%CO2)
    ! deallocate the snow_in
    deallocate(snow_in)
end program TECO