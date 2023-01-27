module mod_ncd_io
    use netcdf
    use mod_data
    implicit none
    
    contains
    subroutine spruce_mip_cmip6Format()
        use netcdf
        implicit none
        ! Daily and monthly
        ! carbon flux (KgC m-2 s-1): gpp, npp, nppLeaf, nppWood, nppRoot, nppOther,
        !              ra, raLeaf, raStem, raRoot, raOther, rMaint, rGrowth, rh
        !              nbp (=gpp - Rh - Ra - other losses)
        !              wetlandCH4, wetlandCH4prod, wetlandCH4cons
        ! carbon pools (KgC m-2): cLeaf, cStem, cRoot, cOther, cLitter (excluding coarse wood debris), cLitterCwd
        !              cSoil, cSoilLevels, cSoilPools (soil organic carbon for each pool), cCH4 (Methane concentration)
        ! Nitrogen flux (KgN m-2 s-1) : fBNF(biological nitrogen fixation), fN2O, fNloss, fNnetmin, fNdep
        ! Nitrogen pools (KgN m-2): nleaf, nStem, nRoot, nOther, nLitter, nLitterCwd, nSoil, nMineral
        ! Energy Fluxes (W m-2): hfls(sensible heat flux), hfss(Latent heat flux), SWnet (Net Shortwave radiation), LWnet(Net Longwave radiation)
        ! Water Fluxes  (Kg m-2 s-1): ec(canopy evaporation), tran(canopy transpiration), es(soil evaporation), hfsbl (snow sublimation), mrro(total runoff),
        !                 mrros (surface runoff), mrrob(subsurface runoff)
        ! other         : mrso (soil moisture in each soil layer, Kg m-2), tsl(soil temperature in each soil layer, K), tsland(surface temperature, K),
        !                 wtd (Water table depth, m), snd (total snow depth, m), lai(m2 m-2) 
        ! ===================================================================================================================================================
        ! carbon fluxes variables
        ! ----------:-----------:----------:-----------------------
        CHARACTER(len=4) :: str_startyr, str_endyr
        
        write(str_startyr,"(I4)")forcing%year(1)
        write(str_endyr,"(I4)")forcing%year(nforcing)
        

        ! hourly outputs
        ! :----------:-----------:----------:-----------:----------:-----------:----------:-----------
        ! hourly GPP
        write(outfile,*) adjustl(trim(outDir_h)),"/gpp_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_gpp_h,"gpp","kgC m-2 s-1", "gross primary productivity")
        ! hourly NPP
        write(outfile,*) adjustl(trim(outDir_h)),"/npp_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_npp_h,"npp","kgC m-2 s-1", "Total net primary productivity")
        ! hourly leaf NPP
        write(outfile,*) adjustl(trim(outDir_h)),"/nppLeaf_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nppLeaf_h,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues")
        ! Hourly wood NPP
        write(outfile,*) adjustl(trim(outDir_h)),"/nppWood_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nppWood_h,"nppWood","kgC m-2 s-1", "NPP allocated to above ground woody tissues")
        ! Hourly stem NPP
        write(outfile,*) adjustl(trim(outDir_h)),"/nppStem_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nppStem_h,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues")
        ! Hourly root NPP
        write(outfile,*) adjustl(trim(outDir_h)),"/nppRoot_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nppRoot_h,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues")
        ! Hourly other NPP
        write(outfile,*) adjustl(trim(outDir_h)),"/nppOther_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nppOther_h,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)")
        ! Hourly ra 
        write(outfile,*) adjustl(trim(outDir_h)),"/ra_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_ra_h,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration")
        ! Hourly leaf ra
        write(outfile,*) adjustl(trim(outDir_h)),"/raLeaf_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_raLeaf_h,"raLeaf","kgC m-2 s-1", "Ra from leaves")
        ! Hourly stem ra
        write(outfile,*) adjustl(trim(outDir_h)),"/raStem_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_raStem_h,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues")
        ! Hourly raRoot_h
        write(outfile,*) adjustl(trim(outDir_h)),"/raRoot_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_raRoot_h,"raRoot","kgC m-2 s-1", "Ra from fine roots")
        ! Hourly raOther_h
        write(outfile,*) adjustl(trim(outDir_h)),"/raOther_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_raOther_h,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)")
        ! Hourly rMaint_h
        write(outfile,*) adjustl(trim(outDir_h)),"/rMaint_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_rMaint_h,"rMaint","kgC m-2 s-1", "Maintenance respiration")
        ! Hourly rGrowth_h                                             ! maintenance respiration and growth respiration
        write(outfile,*) adjustl(trim(outDir_h)),"/rGrowth_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_rGrowth_h,"rGrowth","kgC m-2 s-1", "Growth respiration")
        ! Hourly rh_h
        write(outfile,*) adjustl(trim(outDir_h)),"/rh_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_rh_h,"rh","kgC m-2 s-1", "Heterotrophic respiration rate")
        ! Hourly nbp_h                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        write(outfile,*) adjustl(trim(outDir_h)),"/nbp_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nbp_h,"nbp","kgC m-2 s-1", "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)")
        ! Hourly wetlandCH4_h
        write(outfile,*) adjustl(trim(outDir_h)),"/wetlandCH4_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_wetlandCH4_h,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4")
        ! Hourly wetlandCH4prod_h
        write(outfile,*) adjustl(trim(outDir_h)),"/wetlandCH4prod_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_wetlandCH4prod_h,"wetlandCH4prod","kgC m-2 s-1", "CH4 production")
        ! Hourly wetlandCH4cons_h                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        write(outfile,*) adjustl(trim(outDir_h)),"/wetlandCH4cons_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_wetlandCH4cons_h,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption")

        ! Carbon Pools  (KgC m-2)
        ! Hourly cLeaf_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cLeaf_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cLeaf_h,"cLeaf","kgC m-2", "Carbon biomass in leaves")
        ! Hourly cStem_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cStem_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cStem_h,"cStem","kgC m-2", "Carbon above ground woody biomass")
        ! Hourly cRoot_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cRoot_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cRoot_h,"cRoot","kgC m-2", "Carbon biomass in roots")
        ! Hourly cOther_h                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        write(outfile,*) adjustl(trim(outDir_h)),"/cOther_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cOther_h,"cOther","kgC m-2", "Carbon biomass in other plant organs (reserves, fruits)")
        ! Hourly cLitter_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cLitter_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cLitter_h,"cLitter","kgC m-2", "Carbon in litter (excluding coarse woody debris)")
        ! Hourly cLitterCwd_h                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        write(outfile,*) adjustl(trim(outDir_h)),"/cLitterCwd_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cLitterCwd_h,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris")
        ! Hourly cSoil_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cSoil_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cSoil_h,"cSoil","kgC m-2", "Total soil organic carbon")

        ! Hourly cSoilLevels_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cSoilLevels_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cSoilLevels_h,"cSoilLevels","kgC m-2", "Depth-specific soil organic carbon")
        
        ! Hourly cSoilFast_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cSoilFast_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cSoilFast_h,"cSoilFast","kgC m-2", "Fast soil organic carbon")
        ! Hourly cSoilSlow_h
        write(outfile,*) adjustl(trim(outDir_h)),"/cSoilSlow_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cSoilSlow_h,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon")
        ! Hourly cSoilPassive_h                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        write(outfile,*) adjustl(trim(outDir_h)),"/cSoilPassive_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cSoilPassive_h,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon")
        ! Hourly cCH4_h                                                        ! methane concentration
        write(outfile,*) adjustl(trim(outDir_h)),"/cCH4_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_cCH4_h,"cCH4","kgC m-2 s-1", "Methane concentration")
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! Hourly fBNF_h
        write(outfile,*) adjustl(trim(outDir_h)),"/fBNF_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_fBNF_h,"fBNF","kgN m-2 s-1", "biological nitrogen fixation")
        ! Hourly fN2O_h
        write(outfile,*) adjustl(trim(outDir_h)),"/fN2O_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_fN2O_h,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O")
        ! Hourly fNloss_h
        write(outfile,*) adjustl(trim(outDir_h)),"/fNloss_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_fNloss_h,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching")
        ! Hourly fNnetmin_h
        write(outfile,*) adjustl(trim(outDir_h)),"/fNnetmin_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_fNnetmin_h,"fNnetmin","kgN m-2 s-1", "net mineralization of N")
        ! Hourly fNdep_h                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        write(outfile,*) adjustl(trim(outDir_h)),"/fNdep_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_fNdep_h,"fNdep","kgN m-2 s-1", "Nitrogen deposition")
        
        ! Hourly  Nitrogen pools (kgN m-2)
        ! Hourly nLeaf_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nLeaf_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nLeaf_h,"nLeaf","kgN m-2", "Nitrogen in leaves")
        ! Hourly nStem_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nStem_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nStem_h,"nStem","kgN m-2", "Nitrogen in stems")
        ! Hourly nRoot_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nRoot_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nRoot_h,"nRoot","kgN m-2", "Nirogen in roots")
        ! Hourly nOther_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nOther_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nOther_h,"nOther","kgN m-2", "nitrogen in other plant organs (reserves, fruits)")
        ! Hourly nLitter_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nLitter_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nLitter_h,"nLitter","kgN m-2", "Nitrogen in litter (excluding coarse woody debris)")
        ! Hourly nLitterCwd_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nLitterCwd_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nLitterCwd_h,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris")
        ! Hourly nSoil_h
        write(outfile,*) adjustl(trim(outDir_h)),"/nSoil_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nSoil_h,"nSoil","kgN m-2", "Nitrogen in soil organic matter")
        ! Hourly nMineral_h                    ! nMineral: Mineral nitrogen pool
        write(outfile,*) adjustl(trim(outDir_h)),"/nMineral_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_nMineral_h,"nMineral","kgN m-2", "Mineral nitrogen pool")

        ! Hourly ! energy fluxes (W m-2)
        ! Hourly hfls_h
        write(outfile,*) adjustl(trim(outDir_h)),"/hfls_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_hfls_h,"hfls","W m-2", "Sensible heat flux")
        ! Hourly hfss_h
        write(outfile,*) adjustl(trim(outDir_h)),"/hfss_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_hfss_h,"hfss","W m-2", "Latent heat flux")
        ! Hourly SWnet_h
        write(outfile,*) adjustl(trim(outDir_h)),"/SWnet_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_SWnet_h,"SWnet","W m-2", "Net shortwave radiation")
        ! Hourly LWnet_h                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        write(outfile,*) adjustl(trim(outDir_h)),"/LWnet_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_LWnet_h,"LWnet","W m-2", "Net longwave radiation")

        ! Hourly ! water fluxes (kg m-2 s-1)
        ! Hourly ec_h
        write(outfile,*) adjustl(trim(outDir_h)),"/ec_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_ec_h,"ec","kg m-2 s-1", "Canopy evaporation")
        ! Hourly tran_h
        write(outfile,*) adjustl(trim(outDir_h)),"/tran_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_tran_h,"tran","kg m-2 s-1", "Canopy transpiration")
        ! Hourly es_h                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        write(outfile,*) adjustl(trim(outDir_h)),"/es_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_es_h,"es","kg m-2 s-1", "Soil evaporation")
        ! Hourly hfsbl_h                                                         ! Snow sublimation
        write(outfile,*) adjustl(trim(outDir_h)),"/hfsbl_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_hfsbl_h,"hfsbl","kg m-2 s-1", "Snow sublimation")
        ! Hourly mrro_h
        write(outfile,*) adjustl(trim(outDir_h)),"/mrro_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_mrro_h,"mrro","kg m-2 s-1", "Total runoff")
        ! Hourly mrros_h
        write(outfile,*) adjustl(trim(outDir_h)),"/mrros_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_mrros_h,"mrros","kg m-2 s-1", "Surface runoff")
        ! Hourly mrrob_h                                        ! Total runoff; Surface runoff; Subsurface runoff
        write(outfile,*) adjustl(trim(outDir_h)),"/mrrob_hourly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nHours,all_mrrob_h,"mrrob","kg m-2 s-1", "Subsurface runoff")

        ! daily: 
        ! ---------------------------------------------------------------------
        ! carbon fluxes (KgC m-2 s-1)
        ! Daily gpp_d
        write(outfile,*) adjustl(trim(outDir_d)),"/gpp_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_gpp_d,"gpp","kgC m-2 s-1", "gross primary productivity")
        ! daily NPP
        write(outfile,*) adjustl(trim(outDir_d)),"/npp_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_npp_d,"npp","kgC m-2 s-1", "Total net primary productivity")
        ! daily leaf NPP
        write(outfile,*) adjustl(trim(outDir_d)),"/nppLeaf_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nppLeaf_d,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues")
        ! daily wood NPP
        write(outfile,*) adjustl(trim(outDir_d)),"/nppWood_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nppWood_d,"nppWood","kgC m-2 s-1", "NPP allocated to above ground woody tissues")
        ! daily stem NPP
        write(outfile,*) adjustl(trim(outDir_d)),"/nppStem_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nppStem_d,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues")
        ! daily root NPP
        write(outfile,*) adjustl(trim(outDir_d)),"/nppRoot_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nppRoot_d,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues")
        ! daily other NPP
        write(outfile,*) adjustl(trim(outDir_d)),"/nppOther_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nppOther_d,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)")
        ! daily ra 
        write(outfile,*) adjustl(trim(outDir_d)),"/ra_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_ra_d,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration")
        ! daily leaf ra
        write(outfile,*) adjustl(trim(outDir_d)),"/raLeaf_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_raLeaf_d,"raLeaf","kgC m-2 s-1", "Ra from leaves")
        ! daily stem ra
        write(outfile,*) adjustl(trim(outDir_d)),"/raStem_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_raStem_d,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues")
        ! daily raRoot_d
        write(outfile,*) adjustl(trim(outDir_d)),"/raRoot_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_raRoot_d,"raRoot","kgC m-2 s-1", "Ra from fine roots")
        ! daily raOther_d
        write(outfile,*) adjustl(trim(outDir_d)),"/raOther_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_raOther_d,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)")
        ! daily rMaint_d
        write(outfile,*) adjustl(trim(outDir_d)),"/rMaint_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_rMaint_d,"rMaint","kgC m-2 s-1", "Maintenance respiration")
        ! daily rGrowth_d                                             ! maintenance respiration and growth respiration
        write(outfile,*) adjustl(trim(outDir_d)),"/rGrowth_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_rGrowth_d,"rGrowth","kgC m-2 s-1", "Growth respiration")
        ! daily rh_d
        write(outfile,*) adjustl(trim(outDir_d)),"/rh_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_rh_d,"rh","kgC m-2 s-1", "Heterotrophic respiration rate")
        ! daily nbp_d                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        write(outfile,*) adjustl(trim(outDir_d)),"/nbp_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nbp_d,"nbp","kgC m-2 s-1", "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)")
        ! daily wetlandCH4_d
        write(outfile,*) adjustl(trim(outDir_d)),"/wetlandCH4_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_wetlandCH4_d,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4")
        ! daily wetlandCH4prod_d
        write(outfile,*) adjustl(trim(outDir_d)),"/wetlandCH4prod_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_wetlandCH4prod_d,"wetlandCH4prod","kgC m-2 s-1", "CH4 production")
        ! daily wetlandCH4cons_d                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        write(outfile,*) adjustl(trim(outDir_d)),"/wetlandCH4cons_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_wetlandCH4cons_d,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption")

        ! Carbon Pools  (KgC m-2)
        ! daily cLeaf_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cLeaf_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cLeaf_d,"cLeaf","kgC m-2", "Carbon biomass in leaves")
        ! daily cStem_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cStem_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cStem_d,"cStem","kgC m-2", "Carbon above ground woody biomass")
        ! daily cRoot_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cRoot_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cRoot_d,"cRoot","kgC m-2", "Carbon biomass in roots")
        ! daily cOther_d                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        write(outfile,*) adjustl(trim(outDir_d)),"/cOther_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cOther_d,"cOther","kgC m-2", "Carbon biomass in other plant organs (reserves, fruits)")
        ! daily cLitter_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cLitter_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cLitter_d,"cLitter","kgC m-2", "Carbon in litter (excluding coarse woody debris)")
        ! daily cLitterCwd_d                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        write(outfile,*) adjustl(trim(outDir_d)),"/cLitterCwd_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cLitterCwd_d,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris")
        ! daily cSoil_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cSoil_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cSoil_d,"cSoil","kgC m-2", "Total soil organic carbon")

        ! daily cSoilLevels_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cSoilLevels_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cSoilLevels_d,"cSoilLevels","kgC m-2", "Depth-specific soil organic carbon")
        
        ! daily cSoilFast_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cSoilFast_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cSoilFast_d,"cSoilFast","kgC m-2", "Fast soil organic carbon")
        ! daily cSoilSlow_d
        write(outfile,*) adjustl(trim(outDir_d)),"/cSoilSlow_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cSoilSlow_d,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon")
        ! daily cSoilPassive_d                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        write(outfile,*) adjustl(trim(outDir_d)),"/cSoilPassive_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cSoilPassive_d,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon")
        ! daily cCH4_d                                                        ! methane concentration
        write(outfile,*) adjustl(trim(outDir_d)),"/cCH4_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_cCH4_d,"cCH4","kgC m-2 s-1", "Methane concentration")
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! daily fBNF_d
        write(outfile,*) adjustl(trim(outDir_d)),"/fBNF_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_fBNF_d,"fBNF","kgN m-2 s-1", "biological nitrogen fixation")
        ! daily fN2O_d
        write(outfile,*) adjustl(trim(outDir_d)),"/fN2O_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_fN2O_d,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O")
        ! daily fNloss_d
        write(outfile,*) adjustl(trim(outDir_d)),"/fNloss_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_fNloss_d,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching")
        ! daily fNnetmin_d
        write(outfile,*) adjustl(trim(outDir_d)),"/fNnetmin_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_fNnetmin_d,"fNnetmin","kgN m-2 s-1", "net mineralization of N")
        ! daily fNdep_d                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        write(outfile,*) adjustl(trim(outDir_d)),"/fNdep_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_fNdep_d,"fNdep","kgN m-2 s-1", "Nitrogen deposition")
        
        ! daily  Nitrogen pools (kgN m-2)
        ! daily nLeaf_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nLeaf_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nLeaf_d,"nLeaf","kgN m-2", "Nitrogen in leaves")
        ! daily nStem_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nStem_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nStem_d,"nStem","kgN m-2", "Nitrogen in stems")
        ! daily nRoot_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nRoot_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nRoot_d,"nRoot","kgN m-2", "Nirogen in roots")
        ! daily nOther_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nOther_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nOther_d,"nOther","kgN m-2", "nitrogen in other plant organs (reserves, fruits)")
        ! daily nLitter_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nLitter_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nLitter_d,"nLitter","kgN m-2", "Nitrogen in litter (excluding coarse woody debris)")
        ! daily nLitterCwd_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nLitterCwd_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nLitterCwd_d,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris")
        ! daily nSoil_d
        write(outfile,*) adjustl(trim(outDir_d)),"/nSoil_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nSoil_d,"nSoil","kgN m-2", "Nitrogen in soil organic matter")
        ! daily nMineral_d                    ! nMineral: Mineral nitrogen pool
        write(outfile,*) adjustl(trim(outDir_d)),"/nMineral_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_nMineral_d,"nMineral","kgN m-2", "Mineral nitrogen pool")

        ! daily ! energy fluxes (W m-2)
        ! daily hfls_d
        write(outfile,*) adjustl(trim(outDir_d)),"/hfls_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_hfls_d,"hfls","W m-2", "Sensible heat flux")
        ! daily hfss_d
        write(outfile,*) adjustl(trim(outDir_d)),"/hfss_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_hfss_d,"hfss","W m-2", "Latent heat flux")
        ! daily SWnet_d
        write(outfile,*) adjustl(trim(outDir_d)),"/SWnet_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_SWnet_d,"SWnet","W m-2", "Net shortwave radiation")
        ! daily LWnet_d                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        write(outfile,*) adjustl(trim(outDir_d)),"/LWnet_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_LWnet_d,"LWnet","W m-2", "Net longwave radiation")

        ! daily ! water fluxes (kg m-2 s-1)
        ! daily ec_d
        write(outfile,*) adjustl(trim(outDir_d)),"/ec_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_ec_d,"ec","kg m-2 s-1", "Canopy evaporation")
        ! daily tran_d
        write(outfile,*) adjustl(trim(outDir_d)),"/tran_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_tran_d,"tran","kg m-2 s-1", "Canopy transpiration")
        ! daily es_d                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        write(outfile,*) adjustl(trim(outDir_d)),"/es_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_es_d,"es","kg m-2 s-1", "Soil evaporation")
        ! daily hfsbl_d                                                         ! Snow sublimation
        write(outfile,*) adjustl(trim(outDir_d)),"/hfsbl_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_hfsbl_d,"hfsbl","kg m-2 s-1", "Snow sublimation")
        ! daily mrro_d
        write(outfile,*) adjustl(trim(outDir_d)),"/mrro_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_mrro_d,"mrro","kg m-2 s-1", "Total runoff")
        ! daily mrros_d
        write(outfile,*) adjustl(trim(outDir_d)),"/mrros_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_mrros_d,"mrros","kg m-2 s-1", "Surface runoff")
        ! daily mrrob_d                                        ! Total runoff; Surface runoff; Subsurface runoff
        write(outfile,*) adjustl(trim(outDir_d)),"/mrrob_daily_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nDays,all_mrrob_d,"mrrob","kg m-2 s-1", "Subsurface runoff")

        ! monthly
        ! ---------------------------------------------------------------------
        ! carbon fluxes (KgC m-2 s-1)
        ! monthly gpp_m
        write(outfile,*) adjustl(trim(outDir_m)),"/gpp_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_gpp_m,"gpp","kgC m-2 s-1", "gross primary productivity")
        ! monthly NPP
        write(outfile,*) adjustl(trim(outDir_m)),"/npp_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_npp_m,"npp","kgC m-2 s-1", "Total net primary productivity")
        ! monthly leaf NPP
        write(outfile,*) adjustl(trim(outDir_m)),"/nppLeaf_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nppLeaf_m,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues")
        ! monthly wood NPP
        write(outfile,*) adjustl(trim(outDir_m)),"/nppWood_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nppWood_m,"nppWood","kgC m-2 s-1", "NPP allocated to above ground woody tissues")
        ! monthly stem NPP
        write(outfile,*) adjustl(trim(outDir_m)),"/nppStem_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nppStem_m,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues")
        ! monthly root NPP
        write(outfile,*) adjustl(trim(outDir_m)),"/nppRoot_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nppRoot_m,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues")
        ! monthly other NPP
        write(outfile,*) adjustl(trim(outDir_m)),"/nppOther_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nppOther_m,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)")
        ! monthly ra 
        write(outfile,*) adjustl(trim(outDir_m)),"/ra_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_ra_m,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration")
        ! monthly leaf ra
        write(outfile,*) adjustl(trim(outDir_m)),"/raLeaf_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_raLeaf_m,"raLeaf","kgC m-2 s-1", "Ra from leaves")
        ! monthly stem ra
        write(outfile,*) adjustl(trim(outDir_m)),"/raStem_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_raStem_m,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues")
        ! monthly raRoot_m
        write(outfile,*) adjustl(trim(outDir_m)),"/raRoot_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_raRoot_m,"raRoot","kgC m-2 s-1", "Ra from fine roots")
        ! monthly raOther_m
        write(outfile,*) adjustl(trim(outDir_m)),"/raOther_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_raOther_m,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)")
        ! monthly rMaint_m
        write(outfile,*) adjustl(trim(outDir_m)),"/rMaint_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_rMaint_m,"rMaint","kgC m-2 s-1", "Maintenance respiration")
        ! monthly rGrowth_m                                             ! maintenance respiration and growth respiration
        write(outfile,*) adjustl(trim(outDir_m)),"/rGrowth_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_rGrowth_m,"rGrowth","kgC m-2 s-1", "Growth respiration")
        ! monthly rh_m
        write(outfile,*) adjustl(trim(outDir_m)),"/rh_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_rh_m,"rh","kgC m-2 s-1", "Heterotrophic respiration rate")
        ! monthly nbp_m                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        write(outfile,*) adjustl(trim(outDir_m)),"/nbp_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nbp_m,"nbp","kgC m-2 s-1", "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)")
        ! monthly wetlandCH4_m
        write(outfile,*) adjustl(trim(outDir_m)),"/wetlandCH4_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_wetlandCH4_m,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4")
        ! monthly wetlandCH4prod_m
        write(outfile,*) adjustl(trim(outDir_m)),"/wetlandCH4prod_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_wetlandCH4prod_m,"wetlandCH4prod","kgC m-2 s-1", "CH4 production")
        ! monthly wetlandCH4cons_m                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        write(outfile,*) adjustl(trim(outDir_m)),"/wetlandCH4cons_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_wetlandCH4cons_m,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption")

        ! Carbon Pools  (KgC m-2)
        ! monthly cLeaf_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cLeaf_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cLeaf_m,"cLeaf","kgC m-2", "Carbon biomass in leaves")
        ! monthly cStem_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cStem_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cStem_m,"cStem","kgC m-2", "Carbon above ground woody biomass")
        ! monthly cRoot_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cRoot_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cRoot_m,"cRoot","kgC m-2", "Carbon biomass in roots")
        ! monthly cOther_m                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        write(outfile,*) adjustl(trim(outDir_m)),"/cOther_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cOther_m,"cOther","kgC m-2", "Carbon biomass in other plant organs (reserves, fruits)")
        ! monthly cLitter_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cLitter_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cLitter_m,"cLitter","kgC m-2", "Carbon in litter (excluding coarse woody debris)")
        ! monthly cLitterCwd_m                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        write(outfile,*) adjustl(trim(outDir_m)),"/cLitterCwd_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cLitterCwd_m,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris")
        ! monthly cSoil_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cSoil_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cSoil_m,"cSoil","kgC m-2", "Total soil organic carbon")

        ! monthly cSoilLevels_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cSoilLevels_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cSoilLevels_m,"cSoilLevels","kgC m-2", "Depth-specific soil organic carbon")
        
        ! monthly cSoilFast_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cSoilFast_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cSoilFast_m,"cSoilFast","kgC m-2", "Fast soil organic carbon")
        ! monthly cSoilSlow_m
        write(outfile,*) adjustl(trim(outDir_m)),"/cSoilSlow_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cSoilSlow_m,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon")
        ! monthly cSoilPassive_m                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        write(outfile,*) adjustl(trim(outDir_m)),"/cSoilPassive_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cSoilPassive_m,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon")
        ! monthly cCH4_m                                                        ! methane concentration
        write(outfile,*) adjustl(trim(outDir_m)),"/cCH4_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_cCH4_m,"cCH4","kgC m-2 s-1", "Methane concentration")
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! monthly fBNF_m
        write(outfile,*) adjustl(trim(outDir_m)),"/fBNF_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_fBNF_m,"fBNF","kgN m-2 s-1", "biological nitrogen fixation")
        ! monthly fN2O_m
        write(outfile,*) adjustl(trim(outDir_m)),"/fN2O_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_fN2O_m,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O")
        ! monthly fNloss_m
        write(outfile,*) adjustl(trim(outDir_m)),"/fNloss_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_fNloss_m,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching")
        ! monthly fNnetmin_m
        write(outfile,*) adjustl(trim(outDir_m)),"/fNnetmin_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_fNnetmin_m,"fNnetmin","kgN m-2 s-1", "net mineralization of N")
        ! monthly fNdep_m                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        write(outfile,*) adjustl(trim(outDir_m)),"/fNdep_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_fNdep_m,"fNdep","kgN m-2 s-1", "Nitrogen deposition")
        
        ! monthly  Nitrogen pools (kgN m-2)
        ! monthly nLeaf_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nLeaf_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nLeaf_m,"nLeaf","kgN m-2", "Nitrogen in leaves")
        ! monthly nStem_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nStem_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nStem_m,"nStem","kgN m-2", "Nitrogen in stems")
        ! monthly nRoot_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nRoot_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nRoot_m,"nRoot","kgN m-2", "Nirogen in roots")
        ! monthly nOther_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nOther_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nOther_m,"nOther","kgN m-2", "nitrogen in other plant organs (reserves, fruits)")
        ! monthly nLitter_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nLitter_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nLitter_m,"nLitter","kgN m-2", "Nitrogen in litter (excluding coarse woody debris)")
        ! monthly nLitterCwd_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nLitterCwd_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nLitterCwd_m,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris")
        ! monthly nSoil_m
        write(outfile,*) adjustl(trim(outDir_m)),"/nSoil_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nSoil_m,"nSoil","kgN m-2", "Nitrogen in soil organic matter")
        ! monthly nMineral_m                    ! nMineral: Mineral nitrogen pool
        write(outfile,*) adjustl(trim(outDir_m)),"/nMineral_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_nMineral_m,"nMineral","kgN m-2", "Mineral nitrogen pool")

        ! monthly ! energy fluxes (W m-2)
        ! monthly hfls_m
        write(outfile,*) adjustl(trim(outDir_m)),"/hfls_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_hfls_m,"hfls","W m-2", "Sensible heat flux")
        ! monthly hfss_m
        write(outfile,*) adjustl(trim(outDir_m)),"/hfss_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_hfss_m,"hfss","W m-2", "Latent heat flux")
        ! monthly SWnet_m
        write(outfile,*) adjustl(trim(outDir_m)),"/SWnet_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_SWnet_m,"SWnet","W m-2", "Net shortwave radiation")
        ! monthly LWnet_m                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        write(outfile,*) adjustl(trim(outDir_m)),"/LWnet_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_LWnet_m,"LWnet","W m-2", "Net longwave radiation")

        ! monthly ! water fluxes (kg m-2 s-1)
        ! monthly ec_m
        write(outfile,*) adjustl(trim(outDir_m)),"/ec_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_ec_m,"ec","kg m-2 s-1", "Canopy evaporation")
        ! monthly tran_m
        write(outfile,*) adjustl(trim(outDir_m)),"/tran_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_tran_m,"tran","kg m-2 s-1", "Canopy transpiration")
        ! monthly es_m                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        write(outfile,*) adjustl(trim(outDir_m)),"/es_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_es_m,"es","kg m-2 s-1", "Soil evaporation")
        ! monthly hfsbl_m                                                         ! Snow sublimation
        write(outfile,*) adjustl(trim(outDir_m)),"/hfsbl_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_hfsbl_m,"hfsbl","kg m-2 s-1", "Snow sublimation")
        ! monthly mrro_m
        write(outfile,*) adjustl(trim(outDir_m)),"/mrro_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_mrro_m,"mrro","kg m-2 s-1", "Total runoff")
        ! monthly mrros_m
        write(outfile,*) adjustl(trim(outDir_m)),"/mrros_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_mrros_m,"mrros","kg m-2 s-1", "Surface runoff")
        ! monthly mrrob_m                                        ! Total runoff; Surface runoff; Subsurface runoff
        write(outfile,*) adjustl(trim(outDir_m)),"/mrrob_monthly_TECO-SPRUCE_",adjustl(trim(experiment)),"_", &
            & adjustl(trim(experiment)),"_",adjustl(trim(str_startyr)),"-",adjustl(trim(str_endyr)),".nc"
        call write_nc(outfile,nMonths,all_mrrob_m,"mrrob","kg m-2 s-1", "Subsurface runoff")
        
    end subroutine spruce_mip_cmip6Format
    ! ----------------------------------------------
    subroutine read_nc()

    end subroutine read_nc

    subroutine write_nc(outfile, lenTime, data, varName, unit, description)
        IMPLICIT NONE
        real(kind=4), Dimension(lenTime), intent(in) :: data
        integer(KIND=4) :: ncid, lT_dimid
        integer(kind=4) :: varid
        integer(kind=4), intent(in) :: lenTime
        CHARACTER(LEN=150), INTENT(IN) :: outfile
        CHARACTER(len=*), intent(in) :: varName, unit,description

        
        !Create the netCDF file.
        CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))
        !Define the dimensions.
        CALL check(nf90_def_dim(ncid, "time", lenTime, lT_dimid))
        !Define data variable
        CALL check(nf90_def_var(ncid, varName, NF90_FLOAT, lT_dimid, varid))
        !Add attributes
        CALL check(nf90_put_att(ncid,varid,"units",unit))
        CALL check(nf90_put_att(ncid,varid,"description",description))
        CALL check(nf90_enddef(ncid)) !End Definitions
        !Write Data
        CALL check(nf90_put_var(ncid, varid, data))
        CALL check(nf90_close(ncid))
    end subroutine write_nc
    
    ! check (ever so slightly modified from www.unidata.ucar.edu)
    subroutine check(istatus)
        ! use netcdf
        implicit none
        integer, intent(in) :: istatus
        if(istatus /= nf90_noerr) then
            write(*,*) trim(adjustl(nf90_strerror(istatus)))
        end if
    end subroutine check
end module mod_ncd_io