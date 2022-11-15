module mod_data
    implicit none
    integer :: do_spinup = 0, do_mcmc = 0
    
    ! contant parameters------------------------------------
    real,dimension(3):: tauL,rhoL,rhoS
    real pi,emleaf,emsoil,Rconst,sigma,cpair,Patm,Trefk
    real H2OLv0,airMa,H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta
    real conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,Edjm
    real Entrpy,gam0,gam1,gam2
    ! added for par in methane module   
    real r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi
    ! end of consts parameters------------------------------

    ! ! read parameter from parameter files-----------------
    character(len=150) parafile,commts
    real lat,longi,wsmax,wsmin                  
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmax0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset,Q10,Rl0,Rs0,Rr0
    ! end of read parameters --------------------------------

    integer testtt
    contains
    subroutine get_params()
        ! read the parameters
        implicit none
        parafile=TRIM(parafile)
        ! open and read input file for getting climate data
        open(10,file=parafile,status='old')
        read(10,11)commts
        read(10,*)lat,longi,wsmax,wsmin
        read(10,11)commts
        read(10,*)LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax    
        read(10,11)commts
        read(10,*)SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
        read(10,11)commts
        read(10,*)a1,Ds0,Vcmax0,extkU,xfang,alpha
        read(10,11)commts
        read(10,*)Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,Tau_Micro,Tau_slowSOM,Tau_Passive
        read(10,11)commts
        read(10,*)gddonset,Q10,Rl0,Rs0,Rr0
        ! *** ..int added for pars in methane module
        read(10,11)commts
        read(10,*)r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi                      !this line is for MCMEME
        ! ***     
    11  format(a132)
        close(10)
    end subroutine get_params

    subroutine initialize()
        pi = 3.1415926
        ! physical constants
        tauL(1) = 0.1                  ! leaf transmittance for vis
        rhoL(1) = 0.1                  ! leaf reflectance for vis
        rhoS(1) = 0.1                  ! soil reflectance for vis
        tauL(2) = 0.425                ! for NIR
        rhoL(2) = 0.425                ! for NIR
        rhoS(2) = 0.3                  ! for NIR - later function of soil water content
        tauL(3) = 0.00                 ! for thermal
        rhoL(3) = 0.00                 ! for thermal
        rhoS(3) = 0.00                 ! for thermal
        emleaf  = 0.96
        emsoil  = 0.94
        Rconst  = 8.314                 ! universal gas constant (J/mol)
        sigma   = 5.67e-8                ! Steffan Boltzman constant (W/m2/K4)
        cpair   = 1010.                  ! heat capapcity of air (J/kg/K)
        Patm    = 101325. !1.e5           ! atmospheric pressure  (Pa)
        Trefk   = 293.2                  !reference temp K for Kc, Ko, Rd
        H2OLv0  = 2.501e6               !latent heat H2O (J/kg)
        AirMa   = 29.e-3                 !mol mass air (kg/mol)
        H2OMw   = 18.e-3                 !mol mass H2O (kg/mol)
        chi     = 0.93                     !gbH/gbw
        Dheat   = 21.5e-6                !molecular diffusivity for heat
        ! plant parameters
        gsw0    = 1.0e-2                !g0 for H2O in BWB model
        ! eJmx0 = Vcmx0*2.7            !@20C Leuning 1996 from Wullschleger (1993)
        theta   = 0.9
        wleaf   = 0.01                   !leaf width (m)
        ! thermodynamic parameters for Kc and Ko (Leuning 1990)
        conKc0  = 302.e-6                !mol mol^-1
        conKo0  = 256.e-3                !mol mol^-1
        Ekc     = 59430.                    !J mol^-1
        Eko     = 36000.                    !J mol^-1
        ! Erd = 53000.                    !J mol^-1
        o2ci    = 210.e-3                   !mol mol^-1
        ! thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
        Eavm    = 116300.               !J/mol  (activation energy)
        Edvm    = 202900.               !J/mol  (deactivation energy)
        Eajm    = 79500.                !J/mol  (activation energy) 
        Edjm    = 201000.               !J/mol  (deactivation energy)
        Entrpy  = 650.                !J/mol/K (entropy term, for Jmax & Vcmax)
        ! parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
        gam0    = 28.0e-6               !mol mol^-1 @ 20C = 36.9 @ 25C
        gam1    = .0509
        gam2    = .0010
    end subroutine initialize
    
    subroutine add()
        do_spinup=do_spinup+1
    end subroutine add
end module mod_data