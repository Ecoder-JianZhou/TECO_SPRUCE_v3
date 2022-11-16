module mod_data
    implicit none
    ! run settings 
    integer :: do_spinup = 0, do_mcmc = 0 ! 1 mean run spinup or MCMC
    character(200) climatefile, parafile, commts, snowdepthfile   ! Jian: have the some row number with forcing
    logical, parameter :: do_snow    = .True.
    logical, parameter :: do_soilphy = .True.
    real, parameter :: Ttreat   = 0.
    real, parameter :: CO2treat = 0.
    integer, parameter :: dtimes = 24   ! hourly simulation

    ! for driver --------------
    integer iforcing, nforcing                              ! for cycle
    integer, parameter :: nterms = 11, max_nlines=150000    ! year doy hour Tair Tsoil RH VPD Rain WS PAR CO2
    real,DIMENSION(:), ALLOCATABLE :: snow_in
    type ForcingType
        INTEGER, dimension (:), allocatable :: year
        INTEGER, dimension (:), allocatable :: doy
        INTEGER, dimension (:), allocatable :: hour
        real,    dimension (:), allocatable :: Tair
        real,    dimension (:), allocatable :: Tsoil
        real,    dimension (:), allocatable :: RH           ! RH seems confused in forcing and soil respiration
        real,    dimension (:), allocatable :: VPD
        real,    dimension (:), allocatable :: Rain
        real,    dimension (:), allocatable :: WS
        real,    dimension (:), allocatable :: PAR
        real,    dimension (:), allocatable :: CO2
    end type ForcingType
    type(ForcingType) :: forcing 

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
    real lat,longi,wsmax,wsmin                  
    real LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax
    real SapR,SapS,SLA,GLmax,GRmax,Gsmax,stom_n
    real a1,Ds0,Vcmax0,extkU,xfang,alpha
    real Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C
    real Tau_Micro,Tau_slowSOM,Tau_Passive
    real gddonset,Q10,Rl0,Rs0,Rr0
    ! end of read parameters --------------------------------

    ! initialize parameters ------------------------------------------------------------
    ! parameters for cycle
    integer iyear, iday, ihour
    real radsol, wind, tair, Dair, TairK, co2ca, par, rain, RH, Tsoil
    ! vegetation states and flux
    real GPP, NPP, NSCmin, fnsc, nsc, NSCmax , Vcmx0
    real flait, raero                                   ! Jian: from vegetable to Tsoil_simu
    real RmLeaf,RmStem,RmRoot                           ! maintanence respiration
    real RgLeaf,RgStem,RgRoot                           ! growth respiration 
    real Rgrowth,Rnitrogen,Rmain,Rauto !respirations
    real StemSap,RootSap
    ! soil processes related variables
    real Rsoilab1, Rsoilab2, QLair, QLleaf, QLsoil, Rsoilab3, Rsoilabs      ! calculate in xlayer and used in soil module, QLsoil and Rsoilab3, Rsoilabs seem be calculated both xlayers and soil T dynamic?
    real rhocp, H2OLv, slope, psyc, Cmolar, fw1, Rsoil, rLAI, Hsoil ! seem both in xlayer and soil T dynamic
    real snow_depth, snow_depth_e, snow_dsim
    real,dimension(10):: thksl,wupl,evapl,wcl,FRLEN   ! wsc is the output from soil water module
    real runoff, wsc(10)

    real QC(8) !  leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
    real QN(8),CN0(8),CN(8),OutN(8),QNplant,QNminer
    real N_uptake,N_leach,N_vol,N_fixation,N_deposit,N_fert

    
    
    integer, parameter :: nlayers=10
    real CH4(nlayers),CH4_V(nlayers),CH4V_d(nlayers) 
    real sftmp,Tsnow,Twater,Tice,ice_tw,water_tw 
    ! variables for canopy model
    real evap,transp,ET,G
    real Esoil,Hcrop,ecstot,Anet,DEPH2O,Acanop
    real zwt_d,melt,dcount,dcount_soil
    real,dimension(10):: Tsoill,ice,liq_water
    real zwt,phi
    
    real tau_L,tau_W,tau_R,tau_Micr,tau_Slow,tau_Pass
    real TauC(8)
    real SLAx
    real GLmx,Gsmx,GRmx
    ! for soil conditions
    real WILTPT,FILDCP,infilt
    real stor_use, storage, accumulation
    real SNvcmax,SNgrowth,SNRauto,SNrs
    real LAI,bmroot,bmstem,bmleaf,bmplant,totlivbiom,ht
    
    real N_miner,alphaN
    real NSN, N_deficit
    real, parameter:: times_storage_use=720.   ! 720 hours, 30 days

    real shcap_snow,condu_snow,albedo_snow,resht,thd_snow_depth,b_bound
    ! .. int from soil thermal module
    real diff_snow,diff_s,condu_b
    real depth_ex 
    real infilt_rate
    real fa,fsub,rho_snow,decay_m   
    real fwsoil,topfws,omega 

    real Sps    ! scaling factors for growth

    ! -------------------------------------
    real GDD5
    integer onset !flag of phenological stage
    integer phenoset

    real fwsoil_yr,omega_yr,topfws_yr,diff_yr,diff_d
    real gpp_yr,NPP_yr,NEE_yr,Rh_yr
    real Rh4_yr,Rh5_yr,Rh6_yr,Rh7_yr,Rh8_yr,Ra_yr
    real R_Ntr_yr
    real GL_yr,GR_yr,GW_yr
    real Pool1,Pool2,Pool3,Pool4,Pool5,Pool6,Pool7,Pool8
    real out1_yr,out2_yr,out3_yr,out4_yr,out5_yr,out6_yr,out7_yr,out8_yr

    real runoff_d,runoff_yr,rain_d,rain_yr
    real evap_yr,transp_yr

    real N_up_yr,N_fix_yr,N_dep_yr,N_leach_yr,N_vol_yr

    real soilt_d_simu(11),soilt_d_obs(7),ice_d_simu(10)
    real simuCH4_d, Pro_sum_d, Oxi_sum_d, Fdifu1_d, Ebu_sum_d, Pla_sum_d
    real gpp_d,NPP_d,gpp_ra,NEP_d, NEE_d
    real evap_d,transp_d, Hcanop_d
    real ta, Ts
    real LE_d, RaL, RaS, RaR, Rh_d
    real N_up_d, N_fix_d, N_dep_d, N_leach_d, N_vol_d
    real PAR_d, VPD_d, RECO_d, RLEAV_d, RWOOD_d, RROOT_d
    real GL_d, GW_d, GR_d, LFALL_d, NUP_d, NVOL_d, NLEACH_d
    real NMIN_d, N_LG_d, N_WG_d, N_RG_d, N_LF_d
    real N_WF_d, N_RF_d, WFALL_d, RFALL_d

    integer i
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
        call consts()

        ! plant related parameters
        NSCmin  = 5. 



        fwsoil=1.0
        topfws=1.0
        omega=1.0
        do i=1,10
            wcl(i)=wsmax/100.
        enddo 
        Storage=32.09           !g C/m2
        nsc=85.35
        QC    = (/450.,380.,250.,119.,300.,322.,38340.,23120./) 
        CN0   = (/50.,350.,60.,40.,300.,10.,20.,12./)
        thksl = (/10.,10.,10.,10.,10.,20.,20.,20.,20.,20./)         ! thickness of every soil layer
        FRLEN = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)  ! ratio of roots in every layer, Oak Ridge FACE: Shuang
        ! update: Shuang methane bog species even more shallowly rooted than the tundra
        ! add initials for methane module Shuang version
        CH4_V = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
        CH4   = (/0.0952,0.1232,0.2128,0.3024,0.352,0.8,0.8,0.86,0.86,0.86/)
        ! *** ..int 
        ! add initials for soil thermal dynamics in Yuanyuanversion
        sftmp =-0.
        Tsnow = -20.
        Twater=0.0
        Tice  =0.0

        G=20.5
        Esoil=0.5*G
        snow_dsim =0.575
        dcount=50.
        dcount_soil=50.

        ice_tw = 0.0   
        Tsoill = (/ -0.09, 0.73, 1.3, 1.95, 2.3, 3., 4., 4.5, 5., 5.98/)  ! JJ MS thksl 10 20 30 40 50 70 90 110 130 150...  
        ice    = (/0.1, 0.0, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/) 
            
        liq_water=(/0.01, 0.056, 0.056, 0.056, 0.056, 0.056, 0.056,0.056,0.056,0.056/)    ! unit m
        zwt=0.0
        water_tw=zwt*0.001    
        N_deposit = 2.34/8760. ! Nitrogen input (gN/h/m2, )
        N_fert=0. !5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE)

        ! the unit of residence time is transformed from yearly to hourly
        tauC=(/tau_L,tau_W,tau_R,tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass/)*8760. 
        SLA=SLAx/10000.         ! Convert unit from cm2/g to m2/g
        GLmax=GLmx/8760.    ! growth rates of plant
        GRmax=GRmx/8760.
        Gsmax=GSmx/8760.
        ! end of setting parameters

        ! ===============================================================
        ! cycle  ! skip the following blocks, read input data only.
        ! ===================================================
        ! Initialize parameters and initial state:
        WILTPT=wsmin/100.0
        FILDCP=wsmax/100.0
        ! define soil for export variables for satisfying usage of canopy submodel first time
        ! wscontent=WILTPT
        infilt=0.
        ! gddonset=320.0          
        stor_use=Storage/times_storage_use
        accumulation=0.0
        SNvcmax=1.0
        LAI=LAIMIN
        bmleaf=QC(1)/0.48
        bmstem=QC(2)/0.48
        bmroot=QC(3)/0.48
        bmplant=bmstem+bmroot+bmleaf

        ! initial values of Nitrogen pools and C/N ratio
        alphaN=0.0    ! the transfer of N before littering

        NSN=6.0
        QNminer= 1.2
        N_deficit=0
        CN=CN0
        QN=QC/CN0
        QNplant  =QN(1) + QN(2) + QN(3)

        ! ***********  int initial values of paras used in soil thermal is added here instead of in the pars file ********   
        ! Parameters for soil physical part Yuanyuan 
        ! shcap_snow=690000.  ! refer to csm 4.0 physical constant ! 2090 J/kg/°K with density of 330 kg/m3
        ! ..int 
        shcap_snow=1000000.  ! tuneice worker better
        condu_snow=0.1
        ! condu_b = 0.2  ! yuanyuan int version value
        condu_b = 0.08  ! yuanyuan soil thermal version value  ... int: this par is not sensitive to CWE
        depth_ex=0.05
        ! shcap_snow=700246.3125  ! refer to csm 4.0 physical constant ! 2090 J/kg/°K with density of 330 kg/m3
        ! condu_snow=0.0202
        ! condu_b = 0.0797
        ! depth_ex=0.0472
        
        diff_s=1.
        diff_snow =1.8    ! .. int diffusivity of snow not sensitive for ice
        ! diff_snow =0.018      !tunesnow
        albedo_snow=0.7
        resht=40.
        thd_snow_depth=4.0
        b_bound=100.
        ! b_bound=0.1     !tuneice  not sensitive for ice
        
        infilt_rate= 0.001
        ! infilt_rate= 0.00
        fa = 1
        fsub=0.1
        ! rho_snow=100.
        rho_snow=80.        !tuneice
        decay_m=2.2192      !aging factor on snow melting
    end subroutine initialize

    subroutine init_year()
        implicit none
        GDD5=0.0
        onset=0
        phenoset=0
        diff_yr=0.0
        gpp_yr=0.0
        R_Ntr_yr=0.
        NPP_yr=0.0
        Rh_yr =0.0
        Rh4_yr=0.0
        Rh5_yr=0.0
        Rh6_yr=0.0
        Rh7_yr=0.0
        Rh8_yr=0.0
        Ra_yr =0.0
        GL_yr=0.0
        GW_yr=0.0
        GR_yr=0.0
        Pool1=0.0
        Pool2=0.0
        Pool3=0.0
        Pool4=0.0
        Pool5=0.0
        Pool6=0.0
        Pool7=0.0
        Pool8=0.0
        out1_yr=0.0
        out2_yr=0.0
        out3_yr=0.0
        out4_yr=0.0
        out5_yr=0.0
        out6_yr=0.0
        out7_yr=0.0
        out8_yr=0.0             
        NEE_yr=0.0
        ! water fluxes
        rain_yr=0.0
        transp_yr=0.0
        evap_yr=0.0
        runoff_yr=0.0
        ! Simu_lit=0.
            
        ! Nitrogen fluxes
        N_up_yr=0
        N_fix_yr=0.
        N_dep_yr=0.
        N_leach_yr=0.
        N_vol_yr=0.
        ! ============================== test variable
        fwsoil_yr=0.
        omega_yr=0.
        topfws_yr=0.
        ! hoy=0
    end subroutine init_year

    subroutine init_day()
        implicit none
        ! for daily initials in methane module              
        ! ********* for daily initials in soil thermal module
              
        soilt_d_simu=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
        ice_d_simu=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
        
        soilt_d_obs=(/0.,0.,0.,0.,0.,0.,0./) 
        zwt_d=0.0
        ! --------------------------------------------------
        simuCH4_d=0.0
        Pro_sum_d=0.0
        Oxi_sum_d=0.0
        Fdifu1_d=0.0
        Ebu_sum_d=0.0
        Pla_sum_d=0.0
        CH4V_d= (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)              
        !*****************************               
        !   *** ..int
        ! THE FIRST PART:  coupled canopy and soil model
        diff_d = 0.0
        gpp_d  = 0.0   ! daily
        gpp_ra = 0.0   ! daily
        NPP_d  = 0.0   ! daily
        NEP_d  = 0.0
        NEE_d  = 0.0
        ! rain_d,transp_d,evap_d
        transp_d=0.0   ! daily
        Hcanop_d=0.0   ! daily
        evap_d  =0.0   ! daily
        ta=0.0         ! daily 
        Ts=0.0         ! daily
        rain_d=0.0     ! daily
        runoff_d=0.0    ! daily
        LE_d=0.0
        RaL=0.0
        RaS=0.0
        RaR=0.0
        Rauto=0.0
        Rh_d=0.0
        N_up_d=0.
        N_fix_d=0.
        N_dep_d=0.
        N_leach_d=0.
        N_vol_d=0.
        PAR_d=0.
        VPD_d=0.0
        RECO_d=0.0
        RLEAV_d=0.0 
        RWOOD_d=0.0
        RROOT_d=0.0
        GL_d   =0.0
        GW_d   =0.0
        GR_d   =0.0
        LFALL_d=0.0
        NUP_d=0.0
        NVOL_d=0.
        NLEACH_d=0.0
        NMIN_d=0.0
        N_LG_d=0.0
        N_WG_d=0.0
        N_RG_d=0.0
        N_LF_d=0.0
        N_WF_d=0.0
        N_RF_d=0.0
        WFALL_d=0.0
        RFALL_d=0.0
    end subroutine init_day

    subroutine consts()
        pi = 3.1415926
        ! physical constants
        tauL(1) = 0.1                   ! leaf transmittance for vis
        rhoL(1) = 0.1                   ! leaf reflectance for vis
        rhoS(1) = 0.1                   ! soil reflectance for vis
        tauL(2) = 0.425                 ! for NIR
        rhoL(2) = 0.425                 ! for NIR
        rhoS(2) = 0.3                   ! for NIR - later function of soil water content
        tauL(3) = 0.00                  ! for thermal
        rhoL(3) = 0.00                  ! for thermal
        rhoS(3) = 0.00                  ! for thermal
        emleaf  = 0.96
        emsoil  = 0.94
        Rconst  = 8.314                 ! universal gas constant (J/mol)
        sigma   = 5.67e-8               ! Steffan Boltzman constant (W/m2/K4)
        cpair   = 1010.                 ! heat capapcity of air (J/kg/K)
        Patm    = 101325. !1.e5         ! atmospheric pressure  (Pa)
        Trefk   = 293.2                 ! reference temp K for Kc, Ko, Rd
        H2OLv0  = 2.501e6               ! latent heat H2O (J/kg)
        AirMa   = 29.e-3                ! mol mass air (kg/mol)
        H2OMw   = 18.e-3                ! mol mass H2O (kg/mol)
        chi     = 0.93                  ! gbH/gbw
        Dheat   = 21.5e-6               ! molecular diffusivity for heat
        ! plant parameters
        gsw0    = 1.0e-2                ! g0 for H2O in BWB model
        ! eJmx0 = Vcmx0*2.7             ! @20C Leuning 1996 from Wullschleger (1993)
        theta   = 0.9
        wleaf   = 0.01                  ! leaf width (m)
        ! thermodynamic parameters for Kc and Ko (Leuning 1990)
        conKc0  = 302.e-6               ! mol mol^-1
        conKo0  = 256.e-3               ! mol mol^-1
        Ekc     = 59430.                ! J mol^-1
        Eko     = 36000.                ! J mol^-1
        ! Erd = 53000.                  ! J mol^-1
        o2ci    = 210.e-3               ! mol mol^-1
        ! thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
        Eavm    = 116300.               ! J/mol  (activation energy)
        Edvm    = 202900.               ! J/mol  (deactivation energy)
        Eajm    = 79500.                ! J/mol  (activation energy) 
        Edjm    = 201000.               ! J/mol  (deactivation energy)
        Entrpy  = 650.                  ! J/mol/K (entropy term, for Jmax & Vcmax)
        ! parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
        gam0    = 28.0e-6               ! mol mol^-1 @ 20C = 36.9 @ 25C
        gam1    = .0509
        gam2    = .0010
    end subroutine consts
    
    subroutine get_forcingdata()
        implicit none
        real temp_forcing(nterms, max_nlines)
        integer STAT, COUNT, n, m
        COUNT = 0
        m     = 1
        OPEN(1,FILE=climatefile,status='old',ACTION='read',IOSTAT=STAT)
        read(1,'(a160)') commts
        DO WHILE (.TRUE.)
            COUNT=COUNT+1
            READ(1,*,IOSTAT=STAT) (temp_forcing(n,COUNT), n=1, nterms)
            IF(STAT .NE. 0) EXIT
        ENDDO
        nforcing = COUNT - 1
        CLOSE(1)
        ! initialize the data type of forcing: Tair Tsoil RH VPD Rain WS PAR CO2
        allocate(forcing%year(nforcing))
        allocate(forcing%doy(nforcing))
        allocate(forcing%hour(nforcing))
        allocate(forcing%Tair(nforcing))
        allocate(forcing%Tsoil(nforcing))
        allocate(forcing%RH(nforcing))
        allocate(forcing%VPD(nforcing))
        allocate(forcing%Rain(nforcing))
        allocate(forcing%WS(nforcing))
        allocate(forcing%PAR(nforcing))
        allocate(forcing%CO2(nforcing))
        ! ----------------------------------------------------------------------
        forcing%year  = int(temp_forcing(1,:))
        forcing%doy   = int(temp_forcing(2,:))
        forcing%hour  = int(temp_forcing(3,:))
        forcing%Tair  = temp_forcing(4,:)
        forcing%Tsoil = temp_forcing(5,:)
        forcing%RH    = temp_forcing(6,:)
        forcing%VPD   = temp_forcing(7,:)
        forcing%Rain  = temp_forcing(8,:)
        forcing%WS    = temp_forcing(9,:)
        forcing%PAR   = temp_forcing(10,:)
        forcing%CO2   = temp_forcing(11,:)
    end subroutine get_forcingdata

    subroutine get_snowdepth()
        implicit none
        real temp_snow_depth(max_nlines)
        integer istat1, m

        ! integer m,n,istat1,lines,yr_length
        real snow_depth_read
        integer year,doy,hour

        open(11,file = snowdepthfile, status ='old',ACTION='read', IOSTAT=istat1)
        read(11,'(a160)') commts ! skip 2 lines of input met data file
        m = 0  ! to record the lines in a file
        do
            m=m+1
            read (11,*,IOSTAT=istat1)year,doy,hour,snow_depth_read
            if(istat1<0)exit
            temp_snow_depth(m)=snow_depth_read     
        enddo
        close(11)    ! close snow depth file
        allocate(snow_in(m-1))
        snow_in = temp_snow_depth(1:m-1)
        return
    end subroutine get_snowdepth
end module mod_data