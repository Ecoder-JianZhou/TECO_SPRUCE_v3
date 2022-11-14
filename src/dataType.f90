module mod_data
    implicit none
    integer :: do_spinup = 0, do_mcmc = 0
    
    ! contant parameters
    real, parameter :: pi = 3.1415926
    ! -------------------------------------------
    ! physical constants
    real,dimension(3) :: tauL=[0.1, 0.425, 0.425]   ! leaf transmittance for vis; for NIR; for thermal
    real,dimension(3) :: rhoL=[0.1, 0.425, 0.00]    ! leaf reflectance for vis;   for NIR; for thermal
    real,dimension(3) :: rhoS=[0.1, 0.3,   0.425]   ! soil reflectance for vis;   for NIR; for thermal
    real :: emleaf = 0.96, emsoil=0.94
    real :: Rconst = 8.314                 ! universal gas constant (J/mol)
    real :: sigma  = 5.67e-8                ! Steffan Boltzman constant (W/m2/K4)
    real :: cpair  = 1010.                  ! heat capapcity of air (J/kg/K)
    real :: Patm   = 101325. !1.e5           ! atmospheric pressure  (Pa)
    real :: Trefk  = 293.2                  !reference temp K for Kc, Ko, Rd
    real :: H2OLv0 = 2.501e6               !latent heat H2O (J/kg)
    real :: AirMa  = 29.e-3                 !mol mass air (kg/mol)
    real :: H2OMw  = 18.e-3                 !mol mass H2O (kg/mol)
    real :: chi    = 0.93                     !gbH/gbw
    real :: Dheat  = 21.5e-6                !molecular diffusivity for heat
    ! ! plant parameters
    real :: gsw0   = 1.0e-2                !g0 for H2O in BWB model
    ! real :: eJmx0 = Vcmx0*2.7            !@20C Leuning 1996 from Wullschleger (1993)
    real :: theta  = 0.9
    real :: wleaf  = 0.01                   !leaf width (m)
    ! ! thermodynamic parameters for Kc and Ko (Leuning 1990)
    real :: conKc0 = 302.e-6             !mol mol^-1
    real :: conKo0 = 256.e-3             !mol mol^-1
    real :: Ekc    = 59430.                 !J mol^-1
    real :: Eko    = 36000.                 !J mol^-1
    ! ! Erd = 53000.                     !J mol^-1
    real :: o2ci   = 210.e-3                !mol mol^-1
    ! ! thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
    real :: Eavm   = 116300.               !J/mol  (activation energy)
    real :: Edvm   = 202900.               !J/mol  (deactivation energy)
    real :: Eajm   = 79500.                !J/mol  (activation energy) 
    real :: Edjm   = 201000.               !J/mol  (deactivation energy)
    real :: Entrpy = 650.                !J/mol/K (entropy term, for Jmax & Vcmax)
    ! ! parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
    real :: gam0   = 28.0e-6               !mol mol^-1 @ 20C = 36.9 @ 25C
    real :: gam1   = .0509
    real :: gam2   = .0010
    ! end of consts parameters-----------------------------------------------------------------------

    ! read parameter from parameter files
    real lat,longi,wsmax,wsmin                  
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmax0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset,Q10,Rl0,Rs0,Rr0
    ! end of read parameters ------------------------------------------------------------------------

    integer testtt
    contains
    subroutine initialize()
        gddonset = 100.
        testtt = 666
    end subroutine initialize
    
    subroutine add()
        do_spinup=do_spinup+1
    end subroutine add
end module mod_data