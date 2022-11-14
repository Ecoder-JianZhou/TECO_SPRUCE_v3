module driver
    ! simulation based on the forcing data
    use mod_data
    implicit none
    
    contains
    subroutine teco_simu()
        write(*,*) do_spinup
    end subroutine teco_simu

end module driver