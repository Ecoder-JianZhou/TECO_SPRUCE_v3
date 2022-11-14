program TECO
    use mod_data
    use mod_spinup
    use mod_mcmc
    use driver
    ! to run TECO simulation, spin-up and data simulation
    implicit none

    ! initializations
    call initialize()
    if (do_spinup .eq. 1) call run_spinup() 
    if (do_mcmc .eq. 1) then
        call read_obs()
        call run_mcmc()
    endif
end program TECO