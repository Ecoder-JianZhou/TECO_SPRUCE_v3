module mod_ncd_io
    use netcdf
    use mod_data
    implicit none
    CHARACTER(len=4) :: str_startyr, str_endyr
    
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
        write(str_startyr,"(I4)")forcing%year(1)
        write(str_endyr,"(I4)")forcing%year(nforcing)
        

        ! hourly outputs
        ! :----------:-----------:----------:-----------:----------:-----------:----------:-----------
        ! hourly GPP
        call write_nc(outDir_h,nHours,all_gpp_h,"gpp","kgC m-2 s-1", "gross primary productivity","hourly",1)
        ! hourly NPP
        call write_nc(outDir_h,nHours,all_npp_h,"npp","kgC m-2 s-1", "Total net primary productivity","hourly",1)
        ! hourly leaf NPP
        call write_nc(outDir_h,nHours,all_nppLeaf_h,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues","hourly",1)
        ! Hourly wood NPP
        call write_nc(outDir_h,nHours,all_nppWood_h,"nppWood","kgC m-2 s-1", &
            & "NPP allocated to above ground woody tissues","hourly",1)
        ! Hourly stem NPP
        call write_nc(outDir_h,nHours,all_nppStem_h,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues","hourly",1)
        ! Hourly root NPP
        call write_nc(outDir_h,nHours,all_nppRoot_h,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues","hourly",1)
        ! Hourly other NPP
        call write_nc(outDir_h,nHours,all_nppOther_h,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)","hourly",1)
        ! Hourly ra 
        call write_nc(outDir_h,nHours,all_ra_h,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration","hourly",1)
        ! Hourly leaf ra
        call write_nc(outDir_h,nHours,all_raLeaf_h,"raLeaf","kgC m-2 s-1", "Ra from leaves","hourly",1)
        ! Hourly stem ra
        call write_nc(outDir_h,nHours,all_raStem_h,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues","hourly",1)
        ! Hourly raRoot_h
        call write_nc(outDir_h,nHours,all_raRoot_h,"raRoot","kgC m-2 s-1", "Ra from fine roots","hourly",1)
        ! Hourly raOther_h
        call write_nc(outDir_h,nHours,all_raOther_h,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)","hourly",1)
        ! Hourly rMaint_h
        call write_nc(outDir_h,nHours,all_rMaint_h,"rMaint","kgC m-2 s-1", "Maintenance respiration","hourly",1)
        ! Hourly rGrowth_h                                             ! maintenance respiration and growth respiration
        call write_nc(outDir_h,nHours,all_rGrowth_h,"rGrowth","kgC m-2 s-1", "Growth respiration","hourly",1)
        ! Hourly rh_h
        call write_nc(outDir_h,nHours,all_rh_h,"rh","kgC m-2 s-1", "Heterotrophic respiration rate","hourly",1)
        ! Hourly nbp_h                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        call write_nc(outDir_h,nHours,all_nbp_h,"nbp","kgC m-2 s-1", &
            &"Net Biome productivity (NBP = GPP - Rh - Ra - other losses)","hourly",1)
        ! Hourly wetlandCH4_h
        call write_nc(outDir_h,nHours,all_wetlandCH4_h,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4","hourly",1)
        ! Hourly wetlandCH4prod_h
        call write_nc(outDir_h,nHours,all_wetlandCH4prod_h,"wetlandCH4prod","kgC m-2 s-1", "CH4 production","hourly",1)
        ! Hourly wetlandCH4cons_h                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        call write_nc(outDir_h,nHours,all_wetlandCH4cons_h,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption","hourly",1)

        ! Carbon Pools  (KgC m-2)
        ! Hourly cLeaf_h
        call write_nc(outDir_h,nHours,all_cLeaf_h,"cLeaf","kgC m-2", "Carbon biomass in leaves","hourly",1)
        ! Hourly cStem_h
        call write_nc(outDir_h,nHours,all_cStem_h,"cStem","kgC m-2", "Carbon above ground woody biomass","hourly",1)
        ! Hourly cRoot_h
        call write_nc(outDir_h,nHours,all_cRoot_h,"cRoot","kgC m-2", "Carbon biomass in roots","hourly",1)
        ! Hourly cOther_h                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        call write_nc(outDir_h,nHours,all_cOther_h,"cOther","kgC m-2", &
            & "Carbon biomass in other plant organs (reserves, fruits)","hourly",1)
        ! Hourly cLitter_h
        call write_nc(outDir_h,nHours,all_cLitter_h,"cLitter","kgC m-2", &
            & "Carbon in litter (excluding coarse woody debris)","hourly",1)
        ! Hourly cLitterCwd_h                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        call write_nc(outDir_h,nHours,all_cLitterCwd_h,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris","hourly",1)
        ! Hourly cSoil_h
        call write_nc(outDir_h,nHours,all_cSoil_h,"cSoil","kgC m-2", "Total soil organic carbon","hourly",1)

        ! Hourly cSoilLevels_h
        call write_nc(outDir_h,nHours,all_cSoilLevels_h,"cSoilLevels","kgC m-2", &
            & "Depth-specific soil organic carbon","hourly",nlayers)
        
        ! Hourly cSoilFast_h
        call write_nc(outDir_h,nHours,all_cSoilFast_h,"cSoilFast","kgC m-2", "Fast soil organic carbon","hourly",1)
        ! Hourly cSoilSlow_h
        call write_nc(outDir_h,nHours,all_cSoilSlow_h,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon","hourly",1)
        ! Hourly cSoilPassive_h                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        call write_nc(outDir_h,nHours,all_cSoilPassive_h,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon","hourly",1)
        ! Hourly cCH4_h                                                        ! methane concentration
        call write_nc(outDir_h,nHours,all_cCH4_h,"cCH4","kgC m-2 s-1", "Methane concentration","hourly",nlayers)
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! Hourly fBNF_h
        call write_nc(outDir_h,nHours,all_fBNF_h,"fBNF","kgN m-2 s-1", "biological nitrogen fixation","hourly",1)
        ! Hourly fN2O_h
        call write_nc(outDir_h,nHours,all_fN2O_h,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O","hourly",1)
        ! Hourly fNloss_h
        call write_nc(outDir_h,nHours,all_fNloss_h,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching","hourly",1)
        ! Hourly fNnetmin_h
        call write_nc(outDir_h,nHours,all_fNnetmin_h,"fNnetmin","kgN m-2 s-1", "net mineralization of N","hourly",1)
        ! Hourly fNdep_h                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        call write_nc(outDir_h,nHours,all_fNdep_h,"fNdep","kgN m-2 s-1", "Nitrogen deposition","hourly",1)
        
        ! Hourly  Nitrogen pools (kgN m-2)
        ! Hourly nLeaf_h
        call write_nc(outDir_h,nHours,all_nLeaf_h,"nLeaf","kgN m-2", "Nitrogen in leaves","hourly",1)
        ! Hourly nStem_h
        call write_nc(outDir_h,nHours,all_nStem_h,"nStem","kgN m-2", "Nitrogen in stems","hourly",1)
        ! Hourly nRoot_h
        call write_nc(outDir_h,nHours,all_nRoot_h,"nRoot","kgN m-2", "Nirogen in roots","hourly",1)
        ! Hourly nOther_h
        call write_nc(outDir_h,nHours,all_nOther_h,"nOther","kgN m-2", &
            & "nitrogen in other plant organs (reserves, fruits)","hourly",1)
        ! Hourly nLitter_h
        call write_nc(outDir_h,nHours,all_nLitter_h,"nLitter","kgN m-2", &
            & "Nitrogen in litter (excluding coarse woody debris)","hourly",1)
        ! Hourly nLitterCwd_h
        call write_nc(outDir_h,nHours,all_nLitterCwd_h,"nLitterCwd","kgN m-2", &
            & "Nitrogen in coarse woody debris","hourly",1)
        ! Hourly nSoil_h
        call write_nc(outDir_h,nHours,all_nSoil_h,"nSoil","kgN m-2", "Nitrogen in soil organic matter","hourly",1)
        ! Hourly nMineral_h                    ! nMineral: Mineral nitrogen pool
        call write_nc(outDir_h,nHours,all_nMineral_h,"nMineral","kgN m-2", "Mineral nitrogen pool","hourly",1)

        ! Hourly ! energy fluxes (W m-2)
        ! Hourly hfls_h
        call write_nc(outDir_h,nHours,all_hfls_h,"hfls","W m-2", "Sensible heat flux","hourly",1)
        ! Hourly hfss_h
        call write_nc(outDir_h,nHours,all_hfss_h,"hfss","W m-2", "Latent heat flux","hourly",1)
        ! Hourly SWnet_h
        call write_nc(outDir_h,nHours,all_SWnet_h,"SWnet","W m-2", "Net shortwave radiation","hourly",1)
        ! Hourly LWnet_h                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        call write_nc(outDir_h,nHours,all_LWnet_h,"LWnet","W m-2", "Net longwave radiation","hourly",1)

        ! Hourly ! water fluxes (kg m-2 s-1)
        ! Hourly ec_h
        call write_nc(outDir_h,nHours,all_ec_h,"ec","kg m-2 s-1", "Canopy evaporation","hourly",1)
        ! Hourly tran_h
        call write_nc(outDir_h,nHours,all_tran_h,"tran","kg m-2 s-1", "Canopy transpiration","hourly",1)
        ! Hourly es_h                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        call write_nc(outDir_h,nHours,all_es_h,"es","kg m-2 s-1", "Soil evaporation","hourly",1)
        ! Hourly hfsbl_h                                                         ! Snow sublimation
        call write_nc(outDir_h,nHours,all_hfsbl_h,"hfsbl","kg m-2 s-1", "Snow sublimation","hourly",1)
        ! Hourly mrro_h
        call write_nc(outDir_h,nHours,all_mrro_h,"mrro","kg m-2 s-1", "Total runoff","hourly",1)
        ! Hourly mrros_h
        call write_nc(outDir_h,nHours,all_mrros_h,"mrros","kg m-2 s-1", "Surface runoff","hourly",1)
        ! Hourly mrrob_h                                        ! Total runoff; Surface runoff; Subsurface runoff
        call write_nc(outDir_h,nHours,all_mrrob_h,"mrrob","kg m-2 s-1", "Subsurface runoff","hourly",1)
        ! Hourly Other
        ! Hourly mrso_h
        call write_nc(outDir_h,nHours,all_mrso_h,"mrso","kg m-2", "soil moisture in each soil layer","hourly",nlayers)      ! Kg m-2, soil moisture in each soil layer
        ! Hourly tsl_h 
        call write_nc(outDir_h,nHours,all_tsl_h,"tsl","K", "soil temperature in each soil layer","hourly",nlayers)
        ! Hourly tsland_h
        call write_nc(outDir_h,nHours,all_tsland_h,"tsland","K", "surface temperature","hourly",1)
        ! Hourly wtd_h
        call write_nc(outDir_h,nHours,all_wtd_h,"wtd","m", "Water table depth","hourly",1)
        ! Hourly snd_h
        call write_nc(outDir_h,nHours,all_snd_h,"snd","m", "Total snow depth","hourly",1)
        ! Hourly lai_h
        call write_nc(outDir_h,nHours,all_lai_h,"lai","m2 m-2", "Leaf area index","hourly",1)



        ! daily: 
        ! ---------------------------------------------------------------------
        ! carbon fluxes (KgC m-2 s-1)
        ! Daily gpp_d
        call write_nc(outDir_d,nDays,all_gpp_d,"gpp","kgC m-2 s-1", "gross primary productivity","daily",1)
        ! daily NPP
        call write_nc(outDir_d,nDays,all_npp_d,"npp","kgC m-2 s-1", "Total net primary productivity","daily",1)
        ! daily leaf NPP
        call write_nc(outDir_d,nDays,all_nppLeaf_d,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues","daily",1)
        ! daily wood NPP
        call write_nc(outDir_d,nDays,all_nppWood_d,"nppWood","kgC m-2 s-1", &
            & "NPP allocated to above ground woody tissues","daily",1)
        ! daily stem NPP
        call write_nc(outDir_d,nDays,all_nppStem_d,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues","daily",1)
        ! daily root NPP
        call write_nc(outDir_d,nDays,all_nppRoot_d,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues","daily",1)
        ! daily other NPP
        call write_nc(outDir_d,nDays,all_nppOther_d,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)","daily",1)
        ! daily ra 
        call write_nc(outDir_d,nDays,all_ra_d,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration","daily",1)
        ! daily leaf ra
        call write_nc(outDir_d,nDays,all_raLeaf_d,"raLeaf","kgC m-2 s-1", "Ra from leaves","daily",1)
        ! daily stem ra
        call write_nc(outDir_d,nDays,all_raStem_d,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues","daily",1)
        ! daily raRoot_d
        call write_nc(outDir_d,nDays,all_raRoot_d,"raRoot","kgC m-2 s-1", "Ra from fine roots","daily",1)
        ! daily raOther_d
        call write_nc(outDir_d,nDays,all_raOther_d,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)","daily",1)
        ! daily rMaint_d
        call write_nc(outDir_d,nDays,all_rMaint_d,"rMaint","kgC m-2 s-1", "Maintenance respiration","daily",1)
        ! daily rGrowth_d                                             ! maintenance respiration and growth respiration
        call write_nc(outDir_d,nDays,all_rGrowth_d,"rGrowth","kgC m-2 s-1", "Growth respiration","daily",1)
        ! daily rh_d
        call write_nc(outDir_d,nDays,all_rh_d,"rh","kgC m-2 s-1", "Heterotrophic respiration rate","daily",1)
        ! daily nbp_d                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        call write_nc(outDir_d,nDays,all_nbp_d,"nbp","kgC m-2 s-1", &
            & "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)","daily",1)
        ! daily wetlandCH4_d
        call write_nc(outDir_d,nDays,all_wetlandCH4_d,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4","daily",1)
        ! daily wetlandCH4prod_d
        call write_nc(outDir_d,nDays,all_wetlandCH4prod_d,"wetlandCH4prod","kgC m-2 s-1", "CH4 production","daily",1)
        ! daily wetlandCH4cons_d                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        call write_nc(outDir_d,nDays,all_wetlandCH4cons_d,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption","daily",1)

        ! Carbon Pools  (KgC m-2)
        ! daily cLeaf_d
        call write_nc(outDir_d,nDays,all_cLeaf_d,"cLeaf","kgC m-2", "Carbon biomass in leaves","daily",1)
        ! daily cStem_d
        call write_nc(outDir_d,nDays,all_cStem_d,"cStem","kgC m-2", "Carbon above ground woody biomass","daily",1)
        ! daily cRoot_d
        call write_nc(outDir_d,nDays,all_cRoot_d,"cRoot","kgC m-2", "Carbon biomass in roots","daily",1)
        ! daily cOther_d                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        call write_nc(outDir_d,nDays,all_cOther_d,"cOther","kgC m-2", &
            & "Carbon biomass in other plant organs (reserves, fruits)","daily",1)
        ! daily cLitter_d
        call write_nc(outDir_d,nDays,all_cLitter_d,"cLitter","kgC m-2", &
            & "Carbon in litter (excluding coarse woody debris)","daily",1)
        ! daily cLitterCwd_d                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        call write_nc(outDir_d,nDays,all_cLitterCwd_d,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris","daily",1)
        ! daily cSoil_d
        call write_nc(outDir_d,nDays,all_cSoil_d,"cSoil","kgC m-2", "Total soil organic carbon","daily",1)

        ! daily cSoilLevels_d
        call write_nc(outDir_d,nDays,all_cSoilLevels_d,"cSoilLevels","kgC m-2", &
            & "Depth-specific soil organic carbon","daily",nlayers)
        
        ! daily cSoilFast_d
        call write_nc(outDir_d,nDays,all_cSoilFast_d,"cSoilFast","kgC m-2", "Fast soil organic carbon","daily",1)
        ! daily cSoilSlow_d
        call write_nc(outDir_d,nDays,all_cSoilSlow_d,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon","daily",1)
        ! daily cSoilPassive_d                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        call write_nc(outDir_d,nDays,all_cSoilPassive_d,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon","daily",1)
        ! daily cCH4_d                                                        ! methane concentration
        call write_nc(outDir_d,nDays,all_cCH4_d,"cCH4","kgC m-2 s-1", "Methane concentration","daily",nlayers)
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! daily fBNF_d
        call write_nc(outDir_d,nDays,all_fBNF_d,"fBNF","kgN m-2 s-1", "biological nitrogen fixation","daily",1)
        ! daily fN2O_d
        call write_nc(outDir_d,nDays,all_fN2O_d,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O","daily",1)
        ! daily fNloss_d
        call write_nc(outDir_d,nDays,all_fNloss_d,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching","daily",1)
        ! daily fNnetmin_d
        call write_nc(outDir_d,nDays,all_fNnetmin_d,"fNnetmin","kgN m-2 s-1", "net mineralization of N","daily",1)
        ! daily fNdep_d                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        call write_nc(outDir_d,nDays,all_fNdep_d,"fNdep","kgN m-2 s-1", "Nitrogen deposition","daily",1)
        
        ! daily  Nitrogen pools (kgN m-2)
        ! daily nLeaf_d
        call write_nc(outDir_d,nDays,all_nLeaf_d,"nLeaf","kgN m-2", "Nitrogen in leaves","daily",1)
        ! daily nStem_d
        call write_nc(outDir_d,nDays,all_nStem_d,"nStem","kgN m-2", "Nitrogen in stems","daily",1)
        ! daily nRoot_d
        call write_nc(outDir_d,nDays,all_nRoot_d,"nRoot","kgN m-2", "Nirogen in roots","daily",1)
        ! daily nOther_d
        call write_nc(outDir_d,nDays,all_nOther_d,"nOther","kgN m-2", &
            &"nitrogen in other plant organs (reserves, fruits)","daily",1)
        ! daily nLitter_d
        call write_nc(outDir_d,nDays,all_nLitter_d,"nLitter","kgN m-2",&
            & "Nitrogen in litter (excluding coarse woody debris)","daily",1)
        ! daily nLitterCwd_d
        call write_nc(outDir_d,nDays,all_nLitterCwd_d,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris","daily",1)
        ! daily nSoil_d
        call write_nc(outDir_d,nDays,all_nSoil_d,"nSoil","kgN m-2", "Nitrogen in soil organic matter","daily",1)
        ! daily nMineral_d                    ! nMineral: Mineral nitrogen pool
        call write_nc(outDir_d,nDays,all_nMineral_d,"nMineral","kgN m-2", "Mineral nitrogen pool","daily",1)

        ! daily ! energy fluxes (W m-2)
        ! daily hfls_d
        call write_nc(outDir_d,nDays,all_hfls_d,"hfls","W m-2", "Sensible heat flux","daily",1)
        ! daily hfss_d
        call write_nc(outDir_d,nDays,all_hfss_d,"hfss","W m-2", "Latent heat flux","daily",1)
        ! daily SWnet_d
        call write_nc(outDir_d,nDays,all_SWnet_d,"SWnet","W m-2", "Net shortwave radiation","daily",1)
        ! daily LWnet_d                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        call write_nc(outDir_d,nDays,all_LWnet_d,"LWnet","W m-2", "Net longwave radiation","daily",1)

        ! daily ! water fluxes (kg m-2 s-1)
        ! daily ec_d
        call write_nc(outDir_d,nDays,all_ec_d,"ec","kg m-2 s-1", "Canopy evaporation","daily",1)
        ! daily tran_d
        call write_nc(outDir_d,nDays,all_tran_d,"tran","kg m-2 s-1", "Canopy transpiration","daily",1)
        ! daily es_d                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        call write_nc(outDir_d,nDays,all_es_d,"es","kg m-2 s-1", "Soil evaporation","daily",1)
        ! daily hfsbl_d                                                         ! Snow sublimation
        call write_nc(outDir_d,nDays,all_hfsbl_d,"hfsbl","kg m-2 s-1", "Snow sublimation","daily",1)
        ! daily mrro_d
        call write_nc(outDir_d,nDays,all_mrro_d,"mrro","kg m-2 s-1", "Total runoff","daily",1)
        ! daily mrros_d
        call write_nc(outDir_d,nDays,all_mrros_d,"mrros","kg m-2 s-1", "Surface runoff","daily",1)
        ! daily mrrob_d                                        ! Total runoff; Surface runoff; Subsurface runoff
        call write_nc(outDir_d,nDays,all_mrrob_d,"mrrob","kg m-2 s-1", "Subsurface runoff","daily",1)
        ! daily Other
        ! daily mrso_d
        call write_nc(outDir_d,nDays,all_mrso_d,"mrso","kg m-2", "soil moisture in each soil layer","daily",nlayers)      ! Kg m-2, soil moisture in each soil layer
        ! daily tsl_d 
        call write_nc(outDir_d,nDays,all_tsl_d,"tsl","K", "soil temperature in each soil layer","daily",nlayers)
        ! daily tsland_d
        call write_nc(outDir_d,nDays,all_tsland_d,"tsland","K", "surface temperature","daily",1)
        ! daily wtd_d
        call write_nc(outDir_d,nDays,all_wtd_d,"wtd","m", "Water table depth","daily",1)
        ! daily snd_d
        call write_nc(outDir_d,nDays,all_snd_d,"snd","m", "Total snow depth","daily",1)
        ! daily lai_d
        call write_nc(outDir_d,nDays,all_lai_d,"lai","m2 m-2", "Leaf area index","daily",1)

        ! monthly
        ! ---------------------------------------------------------------------
        ! carbon fluxes (KgC m-2 s-1)
        ! monthly gpp_m
        call write_nc(outDir_m,nMonths,all_gpp_m,"gpp","kgC m-2 s-1", "gross primary productivity","monthly",1)
        ! monthly NPP
        call write_nc(outDir_m,nMonths,all_npp_m,"npp","kgC m-2 s-1", "Total net primary productivity","monthly",1)
        ! monthly leaf NPP
        call write_nc(outDir_m,nMonths,all_nppLeaf_m,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues","monthly",1)
        ! monthly wood NPP
        call write_nc(outDir_m,nMonths,all_nppWood_m,"nppWood","kgC m-2 s-1", &
            & "NPP allocated to above ground woody tissues","monthly",1)
        ! monthly stem NPP
        call write_nc(outDir_m,nMonths,all_nppStem_m,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues","monthly",1)
        ! monthly root NPP
        call write_nc(outDir_m,nMonths,all_nppRoot_m,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues","monthly",1)
        ! monthly other NPP
        call write_nc(outDir_m,nMonths,all_nppOther_m,"nppOther","kgC m-2 s-1", &
            & "NPP allocated to other plant organs (reserves, fruits, exudates)","monthly",1)
        ! monthly ra 
        call write_nc(outDir_m,nMonths,all_ra_m,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration","monthly",1)
        ! monthly leaf ra
        call write_nc(outDir_m,nMonths,all_raLeaf_m,"raLeaf","kgC m-2 s-1", "Ra from leaves","monthly",1)
        ! monthly stem ra
        call write_nc(outDir_m,nMonths,all_raStem_m,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues","monthly",1)
        ! monthly raRoot_m
        call write_nc(outDir_m,nMonths,all_raRoot_m,"raRoot","kgC m-2 s-1", "Ra from fine roots","monthly",1)
        ! monthly raOther_m
        call write_nc(outDir_m,nMonths,all_raOther_m,"raOther","kgC m-2 s-1", &
            & "Ra from other plant organs (reserves, fruits, exudates)","monthly",1)
        ! monthly rMaint_m
        call write_nc(outDir_m,nMonths,all_rMaint_m,"rMaint","kgC m-2 s-1", "Maintenance respiration","monthly",1)
        ! monthly rGrowth_m                                             ! maintenance respiration and growth respiration
        call write_nc(outDir_m,nMonths,all_rGrowth_m,"rGrowth","kgC m-2 s-1", "Growth respiration","monthly",1)
        ! monthly rh_m
        call write_nc(outDir_m,nMonths,all_rh_m,"rh","kgC m-2 s-1", "Heterotrophic respiration rate","monthly",1)
        ! monthly nbp_m                                                    ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        call write_nc(outDir_m,nMonths,all_nbp_m,"nbp","kgC m-2 s-1",&
            & "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)","monthly",1)
        ! monthly wetlandCH4_m
        call write_nc(outDir_m,nMonths,all_wetlandCH4_m,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4","monthly",1)
        ! monthly wetlandCH4prod_m
        call write_nc(outDir_m,nMonths,all_wetlandCH4prod_m,"wetlandCH4prod","kgC m-2 s-1", "CH4 production","monthly",1)
        ! monthly wetlandCH4cons_m                ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        call write_nc(outDir_m,nMonths,all_wetlandCH4cons_m,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption","monthly",1)

        ! Carbon Pools  (KgC m-2)
        ! monthly cLeaf_m
        call write_nc(outDir_m,nMonths,all_cLeaf_m,"cLeaf","kgC m-2", "Carbon biomass in leaves","monthly",1)
        ! monthly cStem_m
        call write_nc(outDir_m,nMonths,all_cStem_m,"cStem","kgC m-2", "Carbon above ground woody biomass","monthly",1)
        ! monthly cRoot_m
        call write_nc(outDir_m,nMonths,all_cRoot_m,"cRoot","kgC m-2", "Carbon biomass in roots","monthly",1)
        ! monthly cOther_m                             ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        call write_nc(outDir_m,nMonths,all_cOther_m,"cOther","kgC m-2", &
            & "Carbon biomass in other plant organs (reserves, fruits)","monthly",1)
        ! monthly cLitter_m
        call write_nc(outDir_m,nMonths,all_cLitter_m,"cLitter","kgC m-2", &
            & "Carbon in litter (excluding coarse woody debris)","monthly",1)
        ! monthly cLitterCwd_m                                         ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        call write_nc(outDir_m,nMonths,all_cLitterCwd_m,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris","monthly",1)
        ! monthly cSoil_m
        call write_nc(outDir_m,nMonths,all_cSoil_m,"cSoil","kgC m-2", "Total soil organic carbon","monthly",1)

        ! monthly cSoilLevels_m
        call write_nc(outDir_m,nMonths,all_cSoilLevels_m,"cSoilLevels","kgC m-2", &
            & "Depth-specific soil organic carbon","monthly",nlayers)
        
        ! monthly cSoilFast_m
        call write_nc(outDir_m,nMonths,all_cSoilFast_m,"cSoilFast","kgC m-2", "Fast soil organic carbon","monthly",1)
        ! monthly cSoilSlow_m
        call write_nc(outDir_m,nMonths,all_cSoilSlow_m,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon","monthly",1)
        ! monthly cSoilPassive_m                            ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        call write_nc(outDir_m,nMonths,all_cSoilPassive_m,"cSoilPassive","kgC m-2 s-1", "Passive soil organic carbon","monthly",1)
        ! monthly cCH4_m                                                        ! methane concentration
        call write_nc(outDir_m,nMonths,all_cCH4_m,"cCH4","kgC m-2 s-1", "Methane concentration","monthly",nlayers)
        
        ! Nitrogen fluxes (kgN m-2 s-1)
        ! monthly fBNF_m
        call write_nc(outDir_m,nMonths,all_fBNF_m,"fBNF","kgN m-2 s-1", "biological nitrogen fixation","monthly",1)
        ! monthly fN2O_m
        call write_nc(outDir_m,nMonths,all_fN2O_m,"fN2O","kgN m-2 s-1", "loss of nitrogen through emission of N2O","monthly",1)
        ! monthly fNloss_m
        call write_nc(outDir_m,nMonths,all_fNloss_m,"fNloss","kgN m-2 s-1", &
            & "Total loss of nitrogen to the atmosphere and from leaching","monthly",1)
        ! monthly fNnetmin_m
        call write_nc(outDir_m,nMonths,all_fNnetmin_m,"fNnetmin","kgN m-2 s-1", "net mineralization of N","monthly",1)
        ! monthly fNdep_m                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        call write_nc(outDir_m,nMonths,all_fNdep_m,"fNdep","kgN m-2 s-1", "Nitrogen deposition","monthly",1)
        
        ! monthly  Nitrogen pools (kgN m-2)
        ! monthly nLeaf_m
        call write_nc(outDir_m,nMonths,all_nLeaf_m,"nLeaf","kgN m-2", "Nitrogen in leaves","monthly",1)
        ! monthly nStem_m
        call write_nc(outDir_m,nMonths,all_nStem_m,"nStem","kgN m-2", "Nitrogen in stems","monthly",1)
        ! monthly nRoot_m
        call write_nc(outDir_m,nMonths,all_nRoot_m,"nRoot","kgN m-2", "Nirogen in roots","monthly",1)
        ! monthly nOther_m
        call write_nc(outDir_m,nMonths,all_nOther_m,"nOther","kgN m-2", &
            &"nitrogen in other plant organs (reserves, fruits)","monthly",1)
        ! monthly nLitter_m
        call write_nc(outDir_m,nMonths,all_nLitter_m,"nLitter","kgN m-2",&
            & "Nitrogen in litter (excluding coarse woody debris)","monthly",1)
        ! monthly nLitterCwd_m
        call write_nc(outDir_m,nMonths,all_nLitterCwd_m,"nLitterCwd","kgN m-2", "Nitrogen in coarse woody debris","monthly",1)
        ! monthly nSoil_m
        call write_nc(outDir_m,nMonths,all_nSoil_m,"nSoil","kgN m-2", "Nitrogen in soil organic matter","monthly",1)
        ! monthly nMineral_m                    ! nMineral: Mineral nitrogen pool
        call write_nc(outDir_m,nMonths,all_nMineral_m,"nMineral","kgN m-2", "Mineral nitrogen pool","monthly",1)

        ! monthly ! energy fluxes (W m-2)
        ! monthly hfls_m
        call write_nc(outDir_m,nMonths,all_hfls_m,"hfls","W m-2", "Sensible heat flux","monthly",1)
        ! monthly hfss_m
        call write_nc(outDir_m,nMonths,all_hfss_m,"hfss","W m-2", "Latent heat flux","monthly",1)
        ! monthly SWnet_m
        call write_nc(outDir_m,nMonths,all_SWnet_m,"SWnet","W m-2", "Net shortwave radiation","monthly",1)
        ! monthly LWnet_m                               ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        call write_nc(outDir_m,nMonths,all_LWnet_m,"LWnet","W m-2", "Net longwave radiation","monthly",1)

        ! monthly ! water fluxes (kg m-2 s-1)
        ! monthly ec_m
        call write_nc(outDir_m,nMonths,all_ec_m,"ec","kg m-2 s-1", "Canopy evaporation","monthly",1)
        ! monthly tran_m
        call write_nc(outDir_m,nMonths,all_tran_m,"tran","kg m-2 s-1", "Canopy transpiration","monthly",1)
        ! monthly es_m                                              ! Canopy evaporation; Canopy transpiration; Soil evaporation
        call write_nc(outDir_m,nMonths,all_es_m,"es","kg m-2 s-1", "Soil evaporation","monthly",1)
        ! monthly hfsbl_m                                                         ! Snow sublimation
        call write_nc(outDir_m,nMonths,all_hfsbl_m,"hfsbl","kg m-2 s-1", "Snow sublimation","monthly",1)
        ! monthly mrro_m
        call write_nc(outDir_m,nMonths,all_mrro_m,"mrro","kg m-2 s-1", "Total runoff","monthly",1)
        ! monthly mrros_m
        call write_nc(outDir_m,nMonths,all_mrros_m,"mrros","kg m-2 s-1", "Surface runoff","monthly",1)
        ! monthly mrrob_m                                        ! Total runoff; Surface runoff; Subsurface runoff
        call write_nc(outDir_m,nMonths,all_mrrob_m,"mrrob","kg m-2 s-1", "Subsurface runoff","monthly",1)
        ! monthly Other
        ! monthly mrso_m
        call write_nc(outDir_m,nMonths,all_mrso_m,"mrso","kg m-2", "soil moisture in each soil layer","monthly",nlayers)      ! Kg m-2, soil moisture in each soil layer
        ! monthly tsl_m 
        call write_nc(outDir_m,nMonths,all_tsl_m,"tsl","K", "soil temperature in each soil layer","monthly",nlayers)
        ! monthly tsland_m
        call write_nc(outDir_m,nMonths,all_tsland_m,"tsland","K", "surface temperature","monthly",1)
        ! monthly wtd_m
        call write_nc(outDir_m,nMonths,all_wtd_m,"wtd","m", "Water table depth","monthly",1)
        ! monthly snd_m
        call write_nc(outDir_m,nMonths,all_snd_m,"snd","m", "Total snow depth","monthly",1)
        ! monthly lai_m
        call write_nc(outDir_m,nMonths,all_lai_m,"lai","m2 m-2", "Leaf area index","monthly",1)
        
    end subroutine spruce_mip_cmip6Format
    ! ----------------------------------------------
    subroutine read_nc()

    end subroutine read_nc

    subroutine write_nc(outfile, lenTime, data, varName, unit, description, freq, nSoilLayer)
        IMPLICIT NONE
        real(kind=4), Dimension(lenTime), intent(in) :: data
        integer(kind=4) :: nSoilLayer
        integer(KIND=4) :: ncid, lT_dimid, dp_dimid
        integer(kind=4) :: varid
        integer(kind=4), intent(in) :: lenTime
        CHARACTER(LEN=*), INTENT(IN) :: outfile, freq
        CHARACTER(len=*), intent(in) :: varName, unit, description
        character(len=:), allocatable :: nc_fileName
        
        allocate(character(len=200+len(outfile)) :: nc_fileName)
        nc_fileName = adjustl(trim(outfile))//"/"//adjustl(trim(varName))//"_"//freq//"_TECO-SPRUCE_"//&
            & adjustl(trim(experiment))//"_"//adjustl(trim(str_startyr))//"-"//adjustl(trim(str_endyr))//".nc"   
        !Create the netCDF file.
        CALL check(nf90_create(nc_fileName, NF90_CLOBBER, ncid))
        !Define the dimensions.
        CALL check(nf90_def_dim(ncid, "time", lenTime, lT_dimid))
        if (nSoilLayer>1)then
            call check(nf90_def_dim(ncid, "depth", nSoilLayer, dp_dimid))
            CALL check(nf90_def_var(ncid, varName, NF90_FLOAT, (/lT_dimid, dp_dimid/), varid))
        else
            CALL check(nf90_def_var(ncid, varName, NF90_FLOAT, lT_dimid, varid))
        endif
        !Define data variable
        
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