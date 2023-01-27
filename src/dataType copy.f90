module mod_data
    implicit none
    ! run settings 
    integer :: do_spinup = 0, do_mcmc = 0                          ! 1 mean run spinup or MCMC
    
    logical, parameter :: do_snow    = .True.                      ! do soil snow process
    logical, parameter :: do_soilphy = .True.                      ! do soil physics
    logical, parameter :: do_matrix  = .True.                      ! do matrix run
    real,    parameter :: Ttreat     = 0.                          ! Temperature treatment, warming in air and soil temperature
    real,    parameter :: CO2treat   = 0.                          ! CO2 treatmant, up to CO2treat, not add to Ca. CO2
    real,    parameter :: N_fert     = 0.                          ! 5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE)
    integer, parameter :: dtimes     = 24                          ! hourly simulation

    ! define some parameters for in/out
    character(len=1500) outdir, outfile                            ! output folder name, output file name (may be changed according to CMIP6 for SPRUCE-MIP)
    character(200) commts, climatefile, parafile, snowdepthfile    ! input file names
    real Simu_dailyflux14(14,10000)                                ! output variables (Jian: need to modify according to CMIP6 for SPURCE-MIP) 
    real Simu_dailyflux14_2023(28,10000)                           ! Jian: Just used for testing Matrix-MIP

    ! parameters for running the loops
    integer nday4out, idayOfnyear, i_record, record_yr(10000), j   ! some variables for cycles
    ! ! parameters for cycle
    integer iyear,  iday, ihour
    real    radsol, wind, co2ca, par, rain, RH
    real    tair, Dair, eairP, TairK, Tsoil                        ! Jian: not sure different between Dair and eairP. eairP means air water vapour pressure 

    ! parameters for forcing data -------------------------------------
    integer iforcing, nforcing                                     ! for cycle
    integer, parameter :: nterms = 11, max_nlines=150000           ! year doy hour Tair Tsoil RH VPD Rain WS PAR CO2
    real,DIMENSION(:), ALLOCATABLE :: snow_in
    type ForcingType
        INTEGER, dimension (:), allocatable :: year
        INTEGER, dimension (:), allocatable :: doy
        INTEGER, dimension (:), allocatable :: hour
        real,    dimension (:), allocatable :: Tair
        real,    dimension (:), allocatable :: Tsoil
        real,    dimension (:), allocatable :: RH                   ! Jian: RH seems confused in forcing and soil respiration
        real,    dimension (:), allocatable :: VPD
        real,    dimension (:), allocatable :: Rain
        real,    dimension (:), allocatable :: WS
        real,    dimension (:), allocatable :: PAR
        real,    dimension (:), allocatable :: CO2
    end type ForcingType
    type(ForcingType) :: forcing 


    ! site-based parameters: read parameter from parameter files ---------------------------------
    real lat, longi
    real wsmax, wsmin
    real LAIMAX, LAIMIN, SLA
    real rdepth
    real Rootmax, Stemmax
    real SapR, SapS
    real GLmax, GRmax, Gsmax
    real stom_n
    real a1
    real Ds0
    real Vcmax0                                             ! Jian: Vcmax0 and Vcmx0 is same? Vcmax0 is Vcmx0 in consts
    real extkU
    real xfang
    real alpha               
    real Tau_Leaf, Tau_Wood, Tau_Root                       ! turnover rate of plant carbon pools : leaf, wood, root  
    real Tau_F,Tau_C                                        ! turnover rate of litter carbon pools: fine, coarse 
    real Tau_Micro, Tau_slowSOM,Tau_Passive                 ! turnover rate of soil carbon pools  : fast, slow, passive 
    real gddonset
    real Q10
    real Rl0, Rs0, Rr0
    ! added for parameters in methane module   
    real r_me
    real Q10pro
    real kCH4
    real Omax
    real CH4_thre
    real Tveg
    real Tpro_me
    real Toxi
    ! end of read parameters --------------------------------

    ! contant parameters------------------------------------
    ! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
    real :: Sps = 1.                                    ! scaling factors for growth, Jian: using in vegetation module. set it to 1
    integer, parameter :: nlayers=10                                ! how many
    ! Jian: like the consts module in old version TECO
    real :: pi      = 3.1415926
    ! physical constants
    real :: tauL(3) = (/0.1, 0.425, 0.00/)          ! leaf transmittance for vis, for NIR, for thermal
    real :: rhoL(3) = (/0.1, 0.425, 0.00/)          ! leaf reflectance for vis, for NIR, for thermal
    real :: rhoS(3) = (/0.1, 0.3,   0.00/)          ! soil reflectance for vis, for NIR, for thermal
    real :: emleaf  = 0.96
    real :: emsoil  = 0.94
    real :: Rconst  = 8.314                         ! universal gas constant (J/mol)
    real :: sigma   = 5.67e-8                       ! Steffan Boltzman constant (W/m2/K4)
    real :: cpair   = 1010.                         ! heat capapcity of air (J/kg/K)
    real :: Patm    = 101325. !1.e5                 ! atmospheric pressure  (Pa)
    real :: Trefk   = 293.2                         ! reference temp K for Kc, Ko, Rd
    real :: H2OLv0  = 2.501e6                       ! latent heat H2O (J/kg)
    real :: AirMa   = 29.e-3                        ! mol mass air (kg/mol)
    real :: H2OMw   = 18.e-3                        ! mol mass H2O (kg/mol)
    real :: chi     = 0.93                          ! gbH/gbw
    real :: Dheat   = 21.5e-6                       ! molecular diffusivity for heat
    ! plant parameters
    real :: gsw0    = 1.0e-2                        ! g0 for H2O in BWB model
    real    eJmx0   != Vcmax0*2.7                   ! @20C Leuning 1996 from Wullschleger (1993)
    real :: theta   = 0.9
    real :: wleaf   = 0.01                          ! leaf width (m)
    ! thermodynamic parameters for Kc and Ko (Leuning 1990)
    real :: conKc0  = 302.e-6                       ! mol mol^-1
    real :: conKo0  = 256.e-3                       ! mol mol^-1
    real :: Ekc     = 59430.                        ! J mol^-1
    real :: Eko     = 36000.                        ! J mol^-1
    ! Erd = 53000.                                  ! J mol^-1
    real :: o2ci    = 210.e-3                       ! mol mol^-1
    ! thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
    real :: Eavm    = 116300.                       ! J/mol  (activation energy)
    real :: Edvm    = 202900.                       ! J/mol  (deactivation energy)
    real :: Eajm    = 79500.                        ! J/mol  (activation energy) 
    real :: Edjm    = 201000.                       ! J/mol  (deactivation energy)
    real :: Entrpy  = 650.                          ! J/mol/K (entropy term, for Jmax & Vcmax)
    ! parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
    real :: gam0    = 28.0e-6                       ! mol mol^-1 @ 20C = 36.9 @ 25C
    real :: gam1    = .0509
    real :: gam2    = .0010
    ! end of consts parameters -------------------------------------------------------------------

    ! initialize parameters ----------------------------------------------------------------------
    ! phonology related parameters
    integer pheno
    real    GDD5
    integer onset                                       !flag of phenological stage
    integer phenoset
    ! ! vegetation states and flux
    real GPP, NPP, NEE, NEP 
    real NPP_L,NPP_W,NPP_R
    real    Rhetero                                         ! total heterotrophic respiration
    real Vcmx0                      
    real NSCmin, fnsc, nsc, NSCmax                      ! none structural carbon pool
    real store, add 
    real L_fall,alpha_L,alpha_W,alpha_R                 ! allocation ratio to Leaf, stem, and Root                      
    real flait, raero                                   ! Jian: from vegetable to Tsoil_simu
    real RmLeaf,RmStem,RmRoot                           ! maintanence respiration
    real RgLeaf,RgStem,RgRoot                           ! growth respiration 
    real Rgrowth,Rnitrogen,Rmain,Rauto                  ! respirations
    real StemSap,RootSap                     
    ! soil processes related variables
    real Rsoilab1, Rsoilab2                             ! calculate in xlayer and used in soil module, QLsoil and Rsoilab3, Rsoilabs seem be calculated both xlayers and soil T dynamic?
    real QLair, QLleaf, QLsoil, Rsoilab3, Rsoilabs      
    real rhocp, H2OLv, slope, psyc, Cmolar, fw1, Rsoil, rLAI, Hsoil ! seem both in xlayer and soil T dynamic
    real snow_depth, snow_depth_e, snow_dsim
    real,dimension(10):: thksl,wcl,FRLEN   ! wsc is the output from soil water module
    real runoff, wsc(10)
    real Rh_pools(5)
    real sublim                             ! snow sublimation
    ! water or energy?
    real Hcanop
    ! methane
    real ProCH4(nlayers), Pro_sum
    real OxiCH4(nlayers), Oxi_sum       !CH4 oxidation
    real Fdifu(nlayers)
    real Ebu_sum,Pla_sum,simuCH4
    real CH4(nlayers),CH4_V(nlayers),CH4V_d(nlayers) 
    ! transfer C
    real QC(8),OutC(8),testout(11)  !  leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
    real QN(8),CN0(8),CN(8),OutN(8),QNplant,QNminer
    real N_uptake,N_leach,N_vol,N_fixation,N_deposit,N_loss,fNnetmin
    real N_leaf,N_wood,N_root
    real N_LF,N_WF,N_RF
    
    
    
    
    real sftmp,Tsnow,Twater,Tice,ice_tw,water_tw 
    ! variables for canopy model
    real evap,transp,ET,G
    real Esoil,Hcrop,ecstot,Anet,DEPH2O,Acanop
    real zwt_d,melt,dcount,dcount_soil
    real,dimension(10):: Tsoill,ice,liq_water
    real zwt,phi
    real Sw
    
    real tau_L,tau_W,tau_R,tau_Micr,tau_Slow,tau_Pass
    real TauC(8)
    real SLAx
    ! real GLmx,Gsmx,GRmx
    ! for soil conditions
    real WILTPT,FILDCP,infilt
    real stor_use, storage, accumulation
    real SNvcmax,SNgrowth,SNRauto,SNrs
    real LAI,bmroot,bmstem,bmleaf,bmplant,totlivbiom,ht
    
    real N_miner,alphaN
    real NSN, N_deficit
    

    real shcap_snow,condu_snow,albedo_snow,resht,thd_snow_depth,b_bound
    ! .. int from soil thermal module
    real diff_snow,diff_s,condu_b
    real depth_ex 
    real infilt_rate
    real fa,fsub,rho_snow,decay_m   
    real fwsoil,topfws,omega 

    

    ! -------------------------------------
    
    ! summary for output or balance calculation variable
    ! hourly: in order to avoid the mistake for hourly or half hour simulation, we define new variable for hourly outputs
    ! -------------------------------------------------------------------------------------------------------------------
    ! carbon fluxes (Kg C m-2 s-1)
    real gpp_h
    real npp_h, nppLeaf_h, nppWood_h, nppStem_h, nppRoot_h, nppOther_h   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
    real ra_h,  raLeaf_h,  raStem_h,  raRoot_h,  raOther_h
    real rMaint_h, rGrowth_h                                             ! maintenance respiration and growth respiration
    real rh_h,  nbp_h                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
    real wetlandCH4_h, wetlandCH4prod_h, wetlandCH4cons_h                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
    ! Carbon Pools  (KgC m-2)
    real cLeaf_h, cStem_h, cRoot_h, cOther_h                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
    real cLitter_h, cLitterCwd_h                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
    real cSoil_h, cSoilLevels_h, cSoilFast_h,cSoilSlow_h,cSoilPassive_h                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
    real cCH4_h(nlayers)                                                          ! methane concentration
    ! Nitrogen fluxes (kgN m-2 s-1)
    real fBNF_h, fN2O_h, fNloss_h, fNnetmin_h, fNdep_h                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
    ! Nitrogen pools (kgN m-2)
    real nLeaf_h, nStem_h, nRoot_h, nOther_h
    real nLitter_h, nLitterCwd_h, nSoil_h, nMineral_h                    ! nMineral: Mineral nitrogen pool
    ! energy fluxes (W m-2)
    real hfls_h, hfss_h, SWnet_h, LWnet_h                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
    ! water fluxes (kg m-2 s-1)
    real ec_h, tran_h, es_h                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
    real hfsbl_h                                                         ! Snow sublimation
    real mrro_h, mrros_h, mrrob_h                                        ! Total runoff; Surface runoff; Subsurface runoff

    ! daily: some daily variables in old TECO need be saved or not? 
    real soilt_d_simu(11),soilt_d_obs(7),ice_d_simu(10)
    real simuCH4_d, Pro_sum_d, Oxi_sum_d, Fdifu1_d, Ebu_sum_d, Pla_sum_d
    real gpp_ra,NEP_d, NEE_d
    real evap_d,transp_d, Hcanop_d
    real ta, Ts, omega_d
    real LE_d, RaL, RaS, RaR
    real N_up_d, N_fix_d, N_dep_d, N_leach_d, N_vol_d
    real PAR_d, VPD_d, RECO_d, RLEAV_d, RWOOD_d, RROOT_d
    real GL_d, GW_d, GR_d, LFALL_d, NUP_d, NVOL_d, NLEACH_d
    real NMIN_d, N_LG_d, N_WG_d, N_RG_d, N_LF_d
    real N_WF_d, N_RF_d, WFALL_d, RFALL_d
    real runoff_d
    real Rsoil_d, ET_d
    ! ---------------------------------------------------------------------
    ! carbon fluxes (Kg C m-2 s-1)
    real gpp_d
    real npp_d, nppLeaf_d, nppWood_d, nppStem_d, nppRoot_d, nppOther_d   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
    real ra_d,  raLeaf_d,  raStem_d,  raRoot_d,  raOther_d
    real rMaint_d, rGrowth_d                                             ! maintenance respiration and growth respiration
    real rh_d,  nbp_d                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
    real wetlandCH4_d, wetlandCH4prod_d, wetlandCH4cons_d                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
    ! Carbon Pools  (KgC m-2)
    real cLeaf_d, cStem_d, cRoot_d, cOther_d                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
    real cLitter_d, cLitterCwd_d                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
    real cSoil_d, cSoilLevels_d, cSoilFast_d,cSoilSlow_d,cSoilPassive_d                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
    real cCH4_d(nlayers)                                                          ! methane concentration
    ! Nitrogen fluxes (kgN m-2 s-1)
    real fBNF_d, fN2O_d, fNloss_d, fNnetmin_d, fNdep_d                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
    ! Nitrogen pools (kgN m-2)
    real nLeaf_d, nStem_d, nRoot_d, nOther_d
    real nLitter_d, nLitterCwd_d, nSoil_d, nMineral_d                    ! nMineral: Mineral nitrogen pool
    ! energy fluxes (W m-2)
    real hfls_d, hfss_d, SWnet_d, LWnet_d                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
    ! water fluxes (kg m-2 s-1)
    real ec_d, tran_d, es_d                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
    real hfsbl_d                                                         ! Snow sublimation
    real mrro_d, mrros_d, mrrob_d                                        ! Total runoff; Surface runoff; Subsurface runoff


    ! monthly
    ! carbon fluxes (Kg C m-2 s-1)
    real gpp_m
    real npp_m, nppLeaf_m, nppWood_m, nppStem_m, nppRoot_m, nppOther_m   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
    real ra_m,  raLeaf_m,  raStem_m,  raRoot_m,  raOther_m
    real rMaint_m, rGrowth_m                                             ! maintenance respiration and growth respiration
    real rh_m,  nbp_m                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
    real wetlandCH4_m, wetlandCH4prod_m, wetlandCH4cons_m                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
    ! Carbon Pools  (KgC m-2)
    real cLeaf_m, cStem_m, cRoot_m, cOther_m                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
    real cLitter_m, cLitterCwd_m                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
    real cSoil_m, cSoilLevels_m, cSoilFast_m,cSoilSlow_m,cSoilPassive_m                           ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
    real cCH4_m(nlayers)                                                          ! methane concentration
    ! Nitrogen fluxes (kgN m-2 s-1)
    real fBNF_m, fN2O_m, fNloss_m, fNnetmin_m, fNdep_m                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
    ! Nitrogen pools (kgN m-2)
    real nLeaf_m, nStem_m, nRoot_m, nOther_m
    real nLitter_m, nLitterCwd_m, nSoil_m, nMineral_m                    ! nMineral: Mineral nitrogen pool
    ! energy fluxes (W m-2)
    real hfls_m, hfss_m, SWnet_m, LWnet_m                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
    ! water fluxes (kg m-2 s-1)
    real ec_m, tran_m, es_m                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
    real hfsbl_m                                                         ! Snow sublimation
    real mrro_m, mrros_m, mrrob_m                                        ! Total runoff; Surface runoff; Subsurface runoff


    ! yearly
    real fwsoil_yr,omega_yr,topfws_yr,diff_yr,diff_d
    real gpp_yr,NPP_yr,NEE_yr,Rh_yr
    real Rh4_yr,Rh5_yr,Rh6_yr,Rh7_yr,Rh8_yr,Ra_yr
    real R_Ntr_yr
    real GL_yr,GR_yr,GW_yr
    real Pool1,Pool2,Pool3,Pool4,Pool5,Pool6,Pool7,Pool8
    real out1_yr,out2_yr,out3_yr,out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
    real runoff_yr,rain_d,rain_yr
    real evap_yr,transp_yr
    real N_up_yr,N_fix_yr,N_dep_yr,N_leach_yr,N_vol_yr
    ! -------------------------------------------------------------------------
    ! carbon fluxes (Kg C m-2 s-1)
    real gpp_y
    real npp_y, nppLeaf_y, nppWood_y, nppStem_y, nppRoot_y, nppOther_y   ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
    real ra_y,  raLeaf_y,  raStem_y,  raRoot_y,  raOther_y
    real rMaint_y, rGrowth_y                                             ! maintenance respiration and growth respiration
    real rh_y,  nbp_y                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
    real wetlandCH4_y, wetlandCH4prod_y, wetlandCH4cons_y                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
    ! Carbon Pools  (KgC m-2)
    real cLeaf_y, cStem_y, cRoot_y, cOther_y                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
    real cLitter_y, cLitterCwd_y                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
    real cSoil_y, cSoilLevels_y, cSoilFast_y, cSoilSlow_y, cSoilPassive_y                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
    real cCH4_y(nlayers)                                                          ! methane concentration
    ! Nitrogen fluxes (kgN m-2 s-1)
    real fBNF_y, fN2O_y, fNloss_y, fNnetmin_y, fNdep_y                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
    ! Nitrogen pools (kgN m-2)
    real nLeaf_y, nStem_y, nRoot_y, nOther_y
    real nLitter_y, nLitterCwd_y, nSoil_y, nMineral_y                    ! nMineral: Mineral nitrogen pool
    ! energy fluxes (W m-2)
    real hfls_y, hfss_y, SWnet_y, LWnet_y                                ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
    ! water fluxes (kg m-2 s-1)
    real ec_y, tran_y, es_y                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
    real hfsbl_y                                                         ! Snow sublimation
    real mrro_y, mrros_y, mrrob_y                                        ! Total runoff; Surface runoff; Subsurface runoff

    ! matrix variables
    real mat_B(8,1), mat_A(8,8), mat_e(8,8), mat_k(8,8), mat_x(8,1), mat_Rh, mat_Rh_d ! for matrix

    integer i

    ! ==============================================================================================
    contains
    subroutine initialize()
        nday4out       = 0
        eJmx0          = Vcmax0*2.7                     ! @20C Leuning 1996 from Wullschleger (1993)
        QC             = (/450.,380.,250.,119.,300.,322.,38340.,23120./)    !  leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
        CN0            = (/50.,350.,60.,40.,300.,10.,20.,12./)
        NSCmin         = 5.                                                 ! none structural carbon pool
        Storage        = 32.09                                              ! g C/m2
        nsc            = 85.35
        stor_use       = Storage/720.                                       ! 720 hours, 30 days
        N_deposit      = 2.34/8760. ! Nitrogen input (gN/h/m2, )
        ! the unit of residence time is transformed from yearly to hourly
        ! tauC           = (/tau_L,tau_W,tau_R,tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass/)*8760. 
        TauC           = (/Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,Tau_Micro,Tau_slowSOM,Tau_Passive/)*8760.
        SLA            = SLAx/10000.          ! Convert unit from cm2/g to m2/g
        GLmax          = GLmax/8760.          ! growth rates of plant. Jian: per year to per hour ?
        GRmax          = GRmax/8760.
        Gsmax          = Gsmax/8760.
        accumulation   = 0.0               ! accumulative storage C?
        SNvcmax        = 1.0
        LAI            = LAIMIN
        bmleaf         = QC(1)/0.48
        bmstem         = QC(2)/0.48
        bmroot         = QC(3)/0.48
        bmplant        = bmstem+bmroot+bmleaf
        ! initial values of Nitrogen pools and C/N ratio
        alphaN         = 0.0    ! the transfer of N before littering
        NSN            = 6.0
        QNminer        = 1.2
        N_deficit      = 0
        CN             = CN0
        QN             = QC/CN0
        QNplant        = QN(1) + QN(2) + QN(3)
        ! -----------------------------------------------------------------------------------------

        ! for soil conditions, physical processes
        thksl          = (/10.,10.,10.,10.,10.,20.,20.,20.,20.,20./)        ! thickness of every soil layer
        FRLEN          = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)  ! ratio of roots in every layer, Oak Ridge FACE: Shuang
        liq_water      = (/0.01, 0.056, 0.056, 0.056, 0.056, 0.056, 0.056,0.056,0.056,0.056/)    ! unit m
        fwsoil         = 1.0                                                ! update in soilwater module
        topfws         = 1.0
        omega          = 1.0
        do i=1,10
            wcl(i)     = wsmax/100.
        enddo 
        zwt            = 0.0
        water_tw       = zwt*0.001 
        WILTPT         = wsmin/100.0
        FILDCP         = wsmax/100.0
        infilt         = 0.
        ! soil thermal dynamics in Yuanyuanversion
        sftmp          = -0.
        Tsnow          = -20.
        Twater         = 0.0
        Tice           = 0.0
        G              = 20.5
        Esoil          = 0.5*G
        snow_dsim      = 0.575
        dcount         = 50.
        dcount_soil    = 50.
        ice_tw         = 0.0   
        Tsoill         = (/ -0.09, 0.73, 1.3, 1.95, 2.3, 3., 4., 4.5, 5., 5.98/)  ! JJ MS thksl 10 20 30 40 50 70 90 110 130 150...  
        ice            = (/0.1, 0.0, 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/) 
        shcap_snow     = 1000000.  ! tuneice worker better
        condu_snow     = 0.1
        condu_b        = 0.08  ! yuanyuan soil thermal version value  ... int: this par is not sensitive to CWE
        depth_ex       = 0.05
        diff_s         = 1.
        diff_snow      = 1.8    ! .. int diffusivity of snow not sensitive for ice
        albedo_snow    = 0.7
        resht          = 40.
        thd_snow_depth = 4.0
        b_bound        = 100.                            ! b_bound=0.1     !tuneice  not sensitive for ice
        infilt_rate    = 0.001
        fa             = 1
        fsub           = 0.1
        ! rho_snow=100.
        rho_snow       = 80.        !tuneice
        decay_m        = 2.2192      !aging factor on snow melting
        !----------------------------------------------------------------
        
        ! methane module. update: Shuang methane bog species even more shallowly rooted than the tundra. add initials for methane module Shuang version
        CH4_V          = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
        CH4            = (/0.0952,0.1232,0.2128,0.3024,0.352,0.8,0.8,0.86,0.86,0.86/)
        do i = 1,8
            mat_x(i,1) = QC(i)
        end do
    end subroutine initialize

    subroutine init_hourly()
        implicit none
        ! Jian: whether it needs to define?
    end subroutine init_hourly

    subroutine init_day()
        implicit none
        ! for daily initials in methane module              
        ! ********* for daily initials in soil thermal module
        soilt_d_simu = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
        ice_d_simu   = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) 
        soilt_d_obs  = (/0.,0.,0.,0.,0.,0.,0./) 
        zwt_d        = 0.0
        ! --------------------------------------------------
        simuCH4_d    = 0.0
        Pro_sum_d    = 0.0
        Oxi_sum_d    = 0.0
        Fdifu1_d     = 0.0
        Ebu_sum_d    = 0.0
        Pla_sum_d    = 0.0
        CH4V_d       = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)              
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
        transp_d = 0.0   ! daily
        Hcanop_d = 0.0   ! daily
        evap_d   = 0.0   ! daily
        ta       = 0.0         ! daily 
        omega_d  = 0.0
        Ts       = 0.0         ! daily
        rain_d   = 0.0     ! daily
        runoff_d = 0.0    ! daily
        LE_d     = 0.0
        RaL      = 0.0
        RaS      = 0.0
        RaR      = 0.0
        Rauto    = 0.0
        Rh_d     = 0.0
        N_up_d   = 0.
        N_fix_d  = 0.
        N_dep_d  = 0.
        N_leach_d= 0.
        N_vol_d  = 0.
        PAR_d    = 0.
        VPD_d    = 0.0
        RECO_d   = 0.0
        RLEAV_d  = 0.0 
        RWOOD_d  = 0.0
        RROOT_d  = 0.0
        GL_d     = 0.0
        GW_d     = 0.0
        GR_d     = 0.0
        LFALL_d  = 0.0
        NUP_d    = 0.0
        NVOL_d   = 0.
        NLEACH_d = 0.0
        NMIN_d   = 0.0
        N_LG_d   = 0.0
        N_WG_d   = 0.0
        N_RG_d   = 0.0
        N_LF_d   = 0.0
        N_WF_d   = 0.0
        N_RF_d   = 0.0
        WFALL_d  = 0.0
        RFALL_d  = 0.0
        mat_Rh_d = 0.0              ! Jian: for matrix model
        ! initialization of daily results according to SPRUCE-MIP
        ! carbon fluxes
        gpp_d             = 0.0
        npp_d             = 0.0
        nppLeaf_d         = 0.0
        nppWood_d         = 0.0
        nppStem_d         = 0.0           ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        nppRoot_d         = 0.0
        nppOther_d        = 0.0
        ra_d              = 0.0
        raLeaf_d          = 0.0
        raStem_d          = 0.0
        raRoot_d          = 0.0
        raOther_d         = 0.0
        rMaint_d          = 0.0           ! maintenance respiration
        rGrowth_d         = 0.0           ! growth respiration
        rh_d              = 0.0           ! heterotrophic respiration
        nbp_d             = 0.0           ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_d      = 0.0           ! wetland net fluxes of CH4
        wetlandCH4prod_d  = 0.0           ! wetland net fluxes of CH4 production
        wetlandCH4cons_d  = 0.0           ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_d           = 0.0
        cStem_d           = 0.0
        cRoot_d           = 0.0
        cOther_d          = 0.0           ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_d         = 0.0           ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_d      = 0.0           ! cLitterCwd: carbon in coarse woody debris
        cSoil_d           = 0.0           ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_d     = 0.0           ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_h       = 0.0           ! cSoilPools (different pools without depth)
        cSoilSlow_h       = 0.0 
        cSoilPassive_h    = 0.0 
        cCH4_d            = 0.0           ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_d            = 0.0           ! fBNF: biological nitrogen fixation;
        fN2O_d            = 0.0           ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_d          = 0.0           ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_d        = 0.0           ! net mineralizaiton
        fNdep_d           = 0.0           ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_d           = 0.0
        nStem_d           = 0.0
        nRoot_d           = 0.0
        nOther_d          = 0.0
        nLitter_d         = 0.0
        nLitterCwd_d      = 0.0
        nSoil_d           = 0.0
        nMineral_d        = 0.0           ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_d            = 0.0           ! Sensible heat flux;
        hfss_d            = 0.0           ! Latent heat flux;
        SWnet_d           = 0.0           ! Net shortwave radiation;
        LWnet_d           = 0.0           ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_d              = 0.0
        tran_d            = 0.0
        es_d              = 0.0           ! Canopy evaporation; Canopy transpiration; Soil evaporation
        hfsbl_d           = 0.0           ! Snow sublimation
        mrro_d            = 0.0
        mrros_d           = 0.0
        mrrob_d           = 0.0           ! Total runoff; Surface runoff; Subsurface runoff
        ! other

        ! not used in SPRUCE-MIP
    end subroutine init_day

    subroutine init_monthly()
        implicit none
        ! Jian: we need add the monthly outputs according to the SPRUCE-MIP's requirement.
        ! carbon result:
        gpp_m             = 0.0
        npp_m             = 0.0
        nppLeaf_m         = 0.0
        nppWood_m         = 0.0
        nppStem_m         = 0.0           ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        nppRoot_m         = 0.0
        nppOther_m        = 0.0
        ra_m              = 0.0
        raLeaf_m          = 0.0
        raStem_m          = 0.0
        raRoot_m          = 0.0
        raOther_m         = 0.0 
        rMaint_m          = 0.0           ! maintenance respiration
        rGrowth_m         = 0.0           ! growth respiration
        rh_m              = 0.0           ! heterotrophic respiration
        nbp_m             = 0.0           ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_m      = 0.0           ! wetland net fluxes of CH4
        wetlandCH4prod_m  = 0.0           ! wetland net fluxes of CH4 production
        wetlandCH4cons_m  = 0.0           ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_m           = 0.0
        cStem_m           = 0.0
        cRoot_m           = 0.0
        cOther_m          = 0.0  ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_m         = 0.0  ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_m      = 0.0  ! cLitterCwd: carbon in coarse woody debris
        cSoil_m           = 0.0  ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_m     = 0.0  ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_m       = 0.0  ! cSoilPools (different pools without depth)
        cSoilSlow_m       = 0.0  
        cSoilPassive_m    = 0.0
        cCH4_m            = 0.0  ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_m            = 0.0  ! fBNF: biological nitrogen fixation;
        fN2O_m            = 0.0  ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_m          = 0.0  ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_m        = 0.0  ! net mineralizaiton
        fNdep_m           = 0.0  ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_m           = 0.0
        nStem_m           = 0.0
        nRoot_m           = 0.0
        nOther_m          = 0.0
        nLitter_m         = 0.0
        nLitterCwd_m      = 0.0
        nSoil_m           = 0.0
        nMineral_m        = 0.0  ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_m            = 0.0  ! Sensible heat flux;
        hfss_m            = 0.0  ! Latent heat flux;
        SWnet_m           = 0.0  ! Net shortwave radiation;
        LWnet_m           = 0.0  !    Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_m              = 0.0
        tran_m            = 0.0
        es_m              = 0.0  ! Canopy evaporation; Canopy transpiration; Soil evaporation
        hfsbl_m           = 0.0  ! Snow sublimation
        mrro_m            = 0.0
        mrros_m           = 0.0
        mrrob_m           = 0.0  ! Total runoff; Surface runoff; Subsurface runoff
        ! other

        ! not used in SPRUCE-MIP
    end subroutine init_monthly

    subroutine init_year()
        implicit none
        GDD5     = 0.0; onset     = 0;   phenoset=0;  diff_yr=0.0; gpp_yr=0.0
        R_Ntr_yr = 0.;  NPP_yr    = 0.0; Rh_yr =0.0;  Rh4_yr=0.0; Rh5_yr=0.0
        Rh6_yr   = 0.0; Rh7_yr    = 0.0; Rh8_yr=0.0;  Ra_yr =0.0; GL_yr=0.0
        GW_yr    = 0.0; GR_yr     = 0.0; Pool1=0.0;   Pool2=0.0; Pool3=0.0
        Pool4    = 0.0; Pool5     = 0.0; Pool6=0.0;   Pool7=0.0; Pool8=0.0
        out1_yr  = 0.0; out2_yr   = 0.0; out3_yr=0.0; out4_yr = 0.0; out5_yr = 0.0
        out6_yr  = 0.0; out7_yr   = 0.0; out8_yr   = 0.0;  NEE_yr    = 0.0
        ! water fluxes
        rain_yr  = 0.0; transp_yr = 0.0; evap_yr   = 0.0; runoff_yr = 0.0
        ! Nitrogen fluxes
        N_up_yr  = 0;   N_fix_yr  = 0.; N_dep_yr=0.; N_leach_yr=0.; N_vol_yr=0.
        ! ============================== test variable
        fwsoil_yr=0.;  omega_yr=0.
        topfws_yr=0.
        ! hoy=0

        ! results according to SPRUCE-MIP
        ! carbon fluxes
        gpp_y             = 0.0
        npp_y             = 0.0
        nppLeaf_y         = 0.0
        nppWood_y         = 0.0   
        nppStem_y         = 0.0 ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        nppRoot_y         = 0.0
        nppOther_y        = 0.0   
        ra_y              = 0.0
        raLeaf_y          = 0.0
        raStem_y          = 0.0
        raRoot_y          = 0.0
        raOther_y         = 0.0
        rMaint_y          = 0.0  ! maintenance respiration
        rGrowth_y         = 0.0  ! growth respiration
        rh_y              = 0.0  ! heterotrophic respiration
        nbp_y             = 0.0  ! NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        wetlandCH4_y      = 0.0  ! wetland net fluxes of CH4
        wetlandCH4prod_y  = 0.0  ! wetland net fluxes of CH4 production
        wetlandCH4cons_y  = 0.0  ! wetland net fluxes of CH4 consumption
        ! Carbon Pools  (KgC m-2)
        cLeaf_y           = 0.0
        cStem_y           = 0.0
        cRoot_y           = 0.0
        cOther_y          = 0.0  ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        cLitter_y         = 0.0  ! litter (excluding coarse woody debris), Jian: fine litter in TECO?
        cLitterCwd_y      = 0.0  ! cLitterCwd: carbon in coarse woody debris
        cSoil_y           = 0.0  ! cSoil: soil organic carbon (Jian: total soil carbon);
        cSoilLevels_y     = 0.0  ! cSoilLevels(depth-specific soil organic carbon, Jian: depth?);
        cSoilFast_y       = 0.0  ! cSoilPools (different pools without depth)
        cSoilSlow_y       = 0.0
        cSoilPassive_y    = 0.0
        cCH4_y            = 0.0  ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        fBNF_y            = 0.0  ! fBNF: biological nitrogen fixation;
        fN2O_y            = 0.0  ! fN2O: loss of nitrogen through emission of N2O;
        fNloss_y          = 0.0  ! fNloss:Total loss of nitrogen to the atmosphere and from leaching;
        fNnetmin_y        = 0.0  ! net mineralizaiton
        fNdep_y           = 0.0  ! deposition of N
        ! Nitrogen pools (kgN m-2)
        nLeaf_y           = 0.0
        nStem_y           = 0.0
        nRoot_y           = 0.0
        nOther_y          = 0.0
        nLitter_y         = 0.0
        nLitterCwd_y      = 0.0
        nSoil_y           = 0.0
        nMineral_y        = 0.0  ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        hfls_y            = 0.0  ! Sensible heat flux;
        hfss_y            = 0.0  ! Latent heat flux;
        SWnet_y           = 0.0  ! Net shortwave radiation;
        LWnet_y           = 0.0  !    Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        ec_y              = 0.0
        tran_y            = 0.0
        es_y              = 0.0  ! Canopy evaporation; Canopy transpiration; Soil evaporation
        hfsbl_y           = 0.0  ! Snow sublimation
        mrro_y            = 0.0
        mrros_y           = 0.0
        mrrob_y           = 0.0  ! Total runoff; Surface runoff; Subsurface runoff
        ! other

        ! not used in SPRUCE-MIP
    end subroutine init_year

    subroutine assign_all_results(hours, days, months, years)
        implicit none
        real, dimension (:), allocatable :: gpp_h

    end subroutine assign_all_results

    ! some functions to get parameters or some input data (eg. forcing data, observation data) 
    !=================================================================
    subroutine get_params()   
        ! Jian: read parameters from the parameter file
        implicit none
        parafile = TRIM(parafile)
        ! open and read input file for getting climate data
        open(10,file=parafile,status='old')
        read(10,11)commts
        read(10,*)lat,longi,wsmax,wsmin
        read(10,11)commts
        read(10,*)LAIMAX,LAIMIN,rdepth,Rootmax,Stemmax    
        read(10,11)commts
        read(10,*)SapR,SapS,SLAx,GLmax,GRmax,Gsmax,stom_n
        read(10,11)commts
        read(10,*)a1,Ds0,Vcmax0,extkU,xfang,alpha
        read(10,11)commts
        read(10,*)Tau_Leaf,Tau_Wood,Tau_Root,Tau_F,Tau_C,Tau_Micro,Tau_slowSOM,Tau_Passive
        read(10,11)commts
        read(10,*)gddonset,Q10,Rl0,Rs0,Rr0
        ! added for pars in methane module
        read(10,11)commts
        read(10,*)r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi                      !this line is for MCMEME     
    11  format(a132)
        close(10)
    end subroutine get_params
    
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