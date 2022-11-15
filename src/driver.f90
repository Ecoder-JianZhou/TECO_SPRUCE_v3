module driver
    ! simulation based on the forcing data
    use mod_data
    use mod_teco
    implicit none
    
    contains
    subroutine teco_simu()
        implicit none
        integer iyear
        do iforcing = 1, nforcing
            if (iforcing .eq. 1) yr0 = forcing_data(1,1)
            iyear = forcing_data(1,iforcing) ! force%year
            if (iyear > yr0) then
                yr0=iyear
                call init_year()
            endif
            iday = forcing_data(2,iforcing)
            StemSap=AMIN1(Stemmax,SapS*bmStem)   ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
            RootSap=AMIN1(Rootmax,SapR*bmRoot)
            NSCmin=5. 
            NSCmax=0.05*(StemSap+RootSap+QC(1))
            if(Ta.gt.5.0)GDD5=GDD5+Ta
            if (do_snow) then 
                if (yr .eq. 1. .and. days .eq. 1.) then
                    ta = -12.85       ! since changed the ta criteria (0. to 1.e-10)) in calculating melt
                    rain_d = 0.        !dbmemo
                endif
                 
                call snow_d(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)                            
                snow_depth_e=snow_dsim
            endif
            call init_day()
            ihour=forcing_data(3,iforcing)
            ! for Duke Forest
            Tair  = input_data(1,iforcing)   ! Tair
            Tsoil = input_data(2,iforcing)    ! SLT
            co2ca = 380.0*1.0E-6
            if (yr .gt. 1)then
                Tair = Tair + Ttreat
                Tsoil = Tsoil + Ttreat
                co2ca=CO2treat*1.0E-6 ! CO2 concentration,ppm-->1.0E-6
            endif                  
            RH   = input_data(3,iforcing)
            Dair = input_data(4,iforcing)       !air water vapour defficit? Unit Pa
            rain  =input_data(5,m)    ! rain fal per hour
            wind  =ABS(input_data(6,m))     ! wind speed m s-1
            PAR   =input_data(7,m)             ! Unit ? umol/s/m-2
            radsol=input_data(7,m)        ! unit ? PAR actually  Jian: or incoming shortwave/longwave radiation?
            co2ca = input_data(8,m)*1.0E-6 ! Jian: add for SPRUCE-MIP 
            ! *** int added for soil thermal/ soil water
            day_mod=mod(m,24)
            if (do_snow) then
                snow_depth=snow_depth_e
            else
                snow_depth=snow_in(m)
            endif
            if (snow_depth .lt. 0.0) snow_depth = 0.0   
                snow_depth = snow_depth*100.   ! change from m to cm  
                ! Rnet=input_data(9,m)
                ! ***
                ! endof Duke Forest

                ! Ajust some unreasonable values
                RH    = AMAX1(0.01,AMIN1(99.99,RH))
                eairP = esat(Tair)*RH/100.             ! Added for SPRUCE, due to lack of VPD data
                Dair  = esat(Tair)-eairP
                radsol=AMAX1(radsol,0.01)
                hoy   =hoy+1
                m=m+1
                if (do_soilphy) then 
                    GOTO 160
                endif
                if(radsol.gt.10.0) then
                    G=-25.0
                else
                    G=20.5											
                endif
                Esoil=0.05*radsol
                if(radsol.LE.10.0) Esoil=0.5*G						
160 continue
                !   ***
                Hcrop=0.1  ! never used in routine
                Ecstot=0.1 ! never used in routine
                Anet=0.1 ! never used in routine
                DepH2O=0.2
                ! for daily mean conditions
                VPD_d=VPD_d+Dair/24./1000.               
                PAR_D=PAR_D+radsol/dtimes                 ! umol photons m-2 s-1   
                ta= ta + tair/24.0             ! sum of a day, for calculating daily mean temperature
                ! write (*,*) 'tair',tair,'ta',ta  !dbmemo
                Ts=Ts+Tsoil/24.0
                rain_d=rain_d+rain
                ! calculating scaling factor of NSC
                if(NSC.le.NSCmin)fnsc=0.0
                if(NSC.ge.NSCmax)fnsc=1.0
                if((NSC.lt.NSCmax).and.(NSC.gt.NSCmin))then 
                    fnsc=(NSC-NSCmin)/(NSCmax-NSCmin)
                endif
                ! update vcmx0 and eJmx0 according to C/N of leaves
                Vcmx0 = Vcmax0*SNvcmax*1.0e-6
                ! eJmx0 = 2.7*Vcmx0  ! original
                eJmx0 = 1.67*Vcmx0 ! Weng 02/21/2011 Medlyn et al. 2002      
                call canopy(gpp,evap,transp,Acanop,Hcanop,Rsoilabs,  & ! outputs
                    &         fwsoil,topfws, &                   ! from soil model
                    &         LAI,Sps,&
                    &         doy,hour,radsol,tair,dair,eairP,    &        ! from climate data file,including 
                    &         wind,rain,&
                    &         Rnet,G,Esoil,Hcrop,Ecstot,Anet,&
                    &         Tsoil,DepH2O,&
                    &         wsmax,wsmin,  &                              !constants specific to soil and plant
                    &         lat,co2ca,a1,Ds0,Vcmx0,extkU,xfang,alpha,&
                    &         stom_n,pi,tauL,rhoL,rhoS,emleaf,emsoil,&
                    &         Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,&
                    &         H2OMw,chi,Dheat,wleaf,gsw0,eJmx0,theta,&
                    &         conKc0,conKo0,Ekc,Eko,o2ci,Eavm,Edvm,Eajm,&
                    &         Edjm,Entrpy,gam0,gam1,gam2,wcl,gddonset,&
                    &         sftmp,Tsoill,testout,ice,&
                    &         water_table_depth,snow_depth,Tsnow,Twater,Tice,water_tw,ice_tw,diff_s,&
                    &         diff_snow,albedo_snow,resht,thd_snow_depth,THKSL,zwt,liq_water,shcap_snow,&
                    &         condu_snow,condu_b,depth_ex,dcount_soil,do_soilphy)

                call soilwater(wsmax,wsmin,rdepth,FRLEN,THKSL,    &   !constants specific to soil/plant
                    &         rain,transp,evap,wcl,runoff,infilt,     &   !inputs
                    &         fwsoil,topfws,omega,wsc,zwt,phi,        &
                    &         liq_water,infilt_rate,melt,ta,day_mod, &
                    &         do_soilphy,snow_depth,ice,testout)                      !outputs

                ET=evap+transp
                rain_yr=rain_yr+rain
                transp_yr=transp_yr+transp
                evap_yr=evap_yr+evap
                runoff_yr=runoff_yr+runoff
                call respiration(LAIMIN,GPP,Tair,Tsoil,DepH2O,&
                    &        Q10,Rl0,Rs0,Rr0,SNRauto,&
                    &        LAI,SLA,bmstem,bmroot,bmleaf,&
                    &        StemSap,RootSap,NSC,fnsc,&
                    &        RmLeaf,RmStem,RmRoot,Rmain)
                ! THE Third Part: update LAI
                call plantgrowth(Tair,omega,GLmax,GRmax,GSmax,&
                    &     LAI,LAIMAX,LAIMIN,SLA,TauC(1),         &    !Tau_L,
                    &     bmleaf,bmroot,bmstem,bmplant,&
                    &     Rootmax,Stemmax,SapS,SapR,&
                    &     StemSap,RootSap,Storage,GDD5,&
                    &     stor_use,onset,accumulation,gddonset,&
                    &     Sps,NSC,fnsc,NSCmin,NSCmax,&
                    &     NSN,CN,CN0,SNgrowth,N_deficit,&
                    &     store,add,L_fall,ht,&
                    &     NPP,alpha_L,alpha_W,alpha_R,&
                    &     RgLeaf,RgStem,RgRoot,Rgrowth)

                ! THE Fourth PART: simulating C influx allocation in pools
                call TCS_CN(Tair,Tsoil,omega,runoff,&
                    &     NPP,alpha_L,alpha_W,alpha_R,L_fall,&
                    &     tauC,QC,OutC,Rh_pools,Rnitrogen,NSC,&
                    &     CNmin,CNmax,NSNmax,NSNmin,alphaN,   &         ! nitrogen
                    &     NSN,N_uptake,N_miner,QN,QNminer,&
                    &     CN,CN0,fnsc,rdepth,N_deficit,&
                    &     N_leaf,N_wood,N_root,N_LF,N_WF,N_RF,&
                    &     N_deposit,N_fixation,N_leach,N_vol,&
                    &     SNvcmax,SNgrowth,SNRauto,SNrs,Q10,&
                    &     tsoill,testout,do_soilphy)
                ! *** ..int 
                ! added tsoil,testout,do_soilphy to TCS_CN
                ! added methane module

                ! write(82,182) zwt, Rh_pools,Tsoil, Ebu_sum_sat, Ebu_sum_unsat
                ! 182     format(5(f15.9,","))       
     
                call methane(Rh_pools,Tsoil,zwt,wsc,     &      !update single value in a hourly loop when MEMCMC=0
                    &     phi,LAIMIN,LAIMAX,           &
                    &     ProCH4,Pro_sum,OxiCH4,Oxi_sum,Fdifu,Ebu_sum,Pla_sum,simuCH4,CH4,CH4_V,   &
                    ! &     ProCH4,Pro_sum,OxiCH4,Oxi_sum,Fdifu,Ebu_sum,Pla_sum,simuCH4,CH4,   &
                    &     r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi,  &
                    &     testout,do_soilphy)       !update single value of Rh_pools,Tsoil,zwt,wsc 
                ! update NSC
                Rauto  =Rmain+Rgrowth+Rnitrogen
                NSC    =NSC+GPP-Rauto-(NPP-add)-store
                Difference = GPP-Rauto-NPP
                if(NSC<0)then
                    bmstem=bmstem+NSC/0.48
                    NPP=NPP+NSC
                    NSN=NSN-NSC/CN(2)
                    NSC=0.
                endif
                GL_d   =GL_d+NPP*alpha_L
                GW_d   =GW_d+NPP*alpha_W
                GR_d   =GR_d+NPP*alpha_R
                LFALL_d=LFALL_d+L_fall
                ! update
                RaLeaf = RgLeaf + RmLeaf
		        RaStem = RgStem + RmStem
		        RaRoot = RgRoot + RmRoot + Rnitrogen
                WFALL_d=WFALL_d+OutC(2) !_wood
                RFALL_d=RFALL_d+OutC(3) !_root
                N_LG_d=N_LG_d+N_leaf
                N_WG_d=N_WG_d+N_wood
                N_RG_d=N_RG_d+N_root
                N_LF_d=N_LF_d+N_LF
                N_WF_d=N_WF_d+N_WF
                N_RF_d=N_RF_d+N_RF

                N_up_d=N_up_d+N_uptake
                N_fix_d=N_fix_d+N_fixation
                N_dep_d=N_dep_d+N_deposit
                N_leach_d=N_leach_d+N_leach
                N_vol_d=N_vol_d+N_vol

                N_up_yr=N_up_yr+N_uptake
                N_fix_yr=N_fix_yr+N_fixation
                N_dep_yr=N_dep_yr+N_deposit
                N_leach_yr=N_leach_yr+N_leach
                N_vol_yr=N_vol_yr+N_vol

                R_Ntr_yr=R_Ntr_yr + Rnitrogen

                ! *** ..int 
                do dlayer=1,10
                    ice_d_simu(dlayer)=ice_d_simu(dlayer)+ice(dlayer) 
                enddo       
                do dlayer=1,11
                    soilt_d_simu(dlayer)=soilt_d_simu(dlayer)+testout(dlayer)  
                    ! first = surface soil temperature 2:11=1:10 layer soil temperatures 
                enddo                    
                do dlayer=1,10
                    CH4V_d(dlayer)=CH4V_d(dlayer)+CH4_V(dlayer) 
                enddo                  
                zwt_d=zwt_d+zwt    ! ..int I doubt it... mean for zwt?     check later  Shuang 
                !   *** 
                ! ==================== test variables
                topfws_yr = topfws_yr+topfws/8760.
                omega_yr=omega_yr+omega/8760.
                fwsoil_yr=fwsoil_yr+fwsoil/8760.

                ! Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
                Rhetero= Rh_pools(1)+Rh_pools(2)+Rh_pools(3) &
                    &    +Rh_pools(4)+Rh_pools(5)
                Rsoil  =Rhetero+RmRoot+RgRoot+Rnitrogen
                NEE=Rauto+Rhetero - GPP
                Q_soil=QC(6) + QC(7) + QC(8)
                bmleaf=QC(1)/0.48
                bmstem=QC(2)/0.48
                bmroot=QC(3)/0.48
                bmplant=bmleaf+bmroot+bmstem
                LAI=bmleaf*SLA
                NMIN_d = NMIN_d+N_miner
                ! output hourly
                Recoh=Rhetero+Rauto
                ETh =ET !*1000.
                Th  =transp !*1000.
                Eh  =evap !*1000.
                INTh=-9999
                VPDh=Dair/1000.
                ROh =runoff !*1000.
                DRAINh=-9999
                LEh =ETh*((2.501-0.00236*Tair)*1000.0)/3600.
                SHh =-9999
                LWh =-9999
                NEP=-NEE
                ! sums of a day
                diff_d=diff_d+difference
                gpp_d=gpp_d + GPP
                gpp_ra=gpp_ra+Rauto
                NPP_d   =NPP_d+NPP
                NEP_d=NEP_d+NEP
                NEE_d=NEE_d+NEE
                RECO_d=RECO_d+Recoh
                Rh_d=  Rh_d + Rhetero
                Ra_d=Reco_d-Rh_d
                RLEAV_d=RLEAV_d+RmLeaf+RgLeaf
                RWOOD_d=RWOOD_d+RmStem+RgStem
                RROOT_d=RROOT_d+RmRoot+RgRoot+Rnitrogen
                Rsoil_d=Rh_d+RROOT_d
                NUP_d=NUP_d+N_uptake
                NVOL_d=NVOL_d+N_vol
                NLEACH_d=NLEACH_d+N_leach
                transp_d=transp_d + transp*(24./dtimes)
                evap_d=evap_d + evap*(24./dtimes)
                ET_d=transp_d + evap_d
                LE_d=LE_d+LEh/24.
                Hcanop_d=Hcanop_d+Hcanop/(24./dtimes)
                runoff_d=runoff_d+runoff
                !   *** .int
                ! added for MEMCMC also for generation of daily methane emission                  
                simuCH4_d=simuCH4_d+simuCH4
                Pro_sum_d=Pro_sum_d+Pro_sum
                Oxi_sum_d=Oxi_sum_d+Oxi_sum
                Fdifu1_d=Fdifu1_d+Fdifu(1)
                Ebu_sum_d=Ebu_sum_d+Ebu_sum
                Pla_sum_d=Pla_sum_d+Pla_sum                                              
                ! ***                  
                ! sum of the whole year
                diff_yr = diff_yr+difference
                gpp_yr=gpp_yr+gpp
                NPP_yr=NPP_yr+NPP
                Rh_yr =Rh_yr +Rhetero
                Ra_yr=Ra_yr+Rauto
                Rh4_yr=Rh4_yr+Rh_pools(1)
                Rh5_yr=Rh5_yr+Rh_pools(2)
                Rh6_yr=Rh6_yr+Rh_pools(3)
                Rh7_yr=Rh7_yr+Rh_pools(4)
                Rh8_yr=Rh8_yr+Rh_pools(5)
                Pool1 = Pool1+QC(1)/8760.
                Pool2 = Pool2+QC(2)/8760.
                Pool3 = Pool3+QC(3)/8760.
                Pool4 = Pool4+QC(4)/8760.
                Pool5 = Pool5+QC(5)/8760.
                Pool6 = Pool6+QC(6)/8760.
                Pool7 = Pool7+QC(7)/8760.
                Pool8 = Pool8+QC(8)/8760.
                out1_yr=out1_yr+OutC(1)
                out2_yr=out2_yr+OutC(2)
                out3_yr=out3_yr+OutC(3)
                out4_yr=out4_yr+OutC(4)
                out5_yr=out5_yr+OutC(5)
                out6_yr=out6_yr+OutC(6)
                out7_yr=out7_yr+OutC(7)
                out8_yr=out8_yr+OutC(8)
                NEE_yr=NEE_yr+NEE
                GL_yr=GL_yr+NPP*alpha_L
                GW_yr=GW_yr+NPP*alpha_W
                GR_yr=GR_yr+NPP*alpha_R
                ! numbering         
                n=n+1
                ! *** .int
                ! added for soil thermal      unknown function check later   Shuang
                if((yr+first_year-1).eq.obs_soilwater(1,k1) .and.    &
                    &     days .eq. obs_soilwater(2,k1) .and.      &
                    &     (i-1).eq. obs_soilwater(3,k1))then
                    Simu_soilwater(1:10,k1)=wcl(1:10)
                    Simu_soiltemp(1:11,k1)=testout
                    Simu_watertable(1,k1)=zwt
                    k1=k1+1
                endif
                ! write(*,*)yr,days,i,gpp,npp
                if(isnan(gpp))then
                    write(*,*)'gpp is nan'
                    return
                endif

            enddo              ! end of dtimes
            if((GDD5.gt.gddonset) .and. phenoset.eq.0) then
                pheno=days
                phenoset=1
            endif

            ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              
            ! output results of canopy and soil models, daily
            ! daily output
            INT_d=-9999
            DRAIN_d=-9999
            SH_d=-9999
            CSOIL_d=QC(6)+QC(7)+QC(8)
            LMA_d=QC(1)/0.48
            NSOIL_d=QN(6)+QN(7)+QN(8)+QNminer

            Rgrowth_d=-9999
            abvLitter=QC(4)  !-9999
            blLitter=QC(5) !-9999
            ! write(*,*)yr,days,LAI,gpp_d,npp_d

            ! if(yr.gt.yrs_eq)then  
            daily=daily+1
  
            Simu_dailyflux(1,daily)=GPP_d ! Leaf
            Simu_dailyflux(2,daily)=NEE_d	! Wood
            Simu_dailyflux(3,daily)=Reco_d	! Coarse roots
            Simu_dailyflux(4,daily)=QC(1)!*1.5
            Simu_dailyflux(5,daily)=GL_yr!*1.5
            Simu_dailyflux(6,daily)=QC(2)!*0.48
            Simu_dailyflux(7,daily)=GW_yr!*0.48
            Simu_dailyflux(8,daily)=QC(3)
            Simu_dailyflux(9,daily)=GR_yr
            Simu_dailyflux(10,daily)=(QC(6)+QC(7)+QC(8))!*13.8   ! Soil
            Simu_dailyflux(11,daily)=pheno   ! Soil
            Simu_dailyflux(12,daily)=LAI !QC(1)/(QC(1)+QC(4)) 

            Simu_dailyflux14(1,daily)=GPP_d 
            Simu_dailyflux14(2,daily)=NEE_d                   
            Simu_dailyflux14(3,daily)=Reco_d	!Rh             
            Simu_dailyflux14(4,daily)=NPP_d!*1.5
            Simu_dailyflux14(5,daily)=Ra_d!*1.5
            Simu_dailyflux14(6,daily)=QC(1)
            Simu_dailyflux14(7,daily)=QC(2)
            Simu_dailyflux14(8,daily)=QC(3)
            Simu_dailyflux14(9,daily)=QC(4)
            Simu_dailyflux14(10,daily)=QC(5)
            Simu_dailyflux14(11,daily)=QC(6)
            Simu_dailyflux14(12,daily)=QC(7)
            Simu_dailyflux14(13,daily)=QC(8)!*0.48
            Simu_dailyflux14(14,daily)=Rh_d
                
            ! Simu_dailywater(1,daily)= wsc(1)        ! not aggregated to daily, value should represents 23:00
            ! Simu_dailywater(2,daily)= wsc(2)
            ! Simu_dailywater(3,daily)= wsc(3)
            ! Simu_dailywater(4,daily)= wsc(4)
            ! Simu_dailywater(5,daily)= wsc(5)
            ! Simu_dailywater(6,daily)= wsc(6)        ! not aggregated to daily, value should represents 23:00
            ! Simu_dailywater(7,daily)= wsc(7)
            ! Simu_dailywater(8,daily)= wsc(8)
            ! Simu_dailywater(9,daily)= wsc(9)
            ! Simu_dailywater(10,daily)= wsc(10)
            !            
            ! if (do_soilphy)
            ! write (*,*) 'simuwcl',wcl(1)           !dbice
            Simu_dailywater(1,daily)= wcl(1)        ! not aggregated to daily, value should represents 23:00
            Simu_dailywater(2,daily)= wcl(2)
            Simu_dailywater(3,daily)= wcl(3)
            Simu_dailywater(4,daily)= wcl(4)
            Simu_dailywater(5,daily)= wcl(5)
            Simu_dailywater(6,daily)= wcl(6)        ! not aggregated to daily, value should represents 23:00
            Simu_dailywater(7,daily)= wcl(7)
            Simu_dailywater(8,daily)= wcl(8)
            Simu_dailywater(9,daily)= wcl(9)
            Simu_dailywater(10,daily)= wcl(10)
            Simu_dailywater(11,daily)= liq_water(1)        ! not aggregated to daily, value should represents 23:00
            Simu_dailywater(12,daily)= liq_water(2)
            Simu_dailywater(13,daily)= liq_water(3)
            Simu_dailywater(14,daily)= liq_water(4)
            Simu_dailywater(15,daily)= liq_water(5)
            Simu_dailywater(16,daily)= liq_water(6)        ! not aggregated to daily, value should represents 23:00
            Simu_dailywater(17,daily)= liq_water(7)
            Simu_dailywater(18,daily)= liq_water(8)
            Simu_dailywater(19,daily)= liq_water(9)
            Simu_dailywater(20,daily)= liq_water(10)
            Simu_dailywater(21,daily)= ice(1)        ! not aggregated to daily, value should represents 23:00
            Simu_dailywater(22,daily)= ice(2)
            Simu_dailywater(23,daily)= ice(3)
            Simu_dailywater(24,daily)= ice(4)
            Simu_dailywater(25,daily)= ice(5)
            Simu_dailywater(26,daily)= ice(6)        ! not aggregated to daily, value should represents 23:00
            Simu_dailywater(27,daily)= ice(7)
            Simu_dailywater(28,daily)= ice(8)
            Simu_dailywater(29,daily)= ice(9)
            Simu_dailywater(30,daily)= ice(10)            
            Simu_dailywater(31,daily)= zwt
            ! *** ..int methane           
            Simu_dailyCH4(1,daily)=simuCH4_d
            Simu_dailyCH4(2,daily)=Pro_sum_d
            Simu_dailyCH4(3,daily)=Oxi_sum_d
            Simu_dailyCH4(4,daily)=Fdifu1_d
            Simu_dailyCH4(5,daily)=Ebu_sum_d
            Simu_dailyCH4(6,daily)=Pla_sum_d
            Simu_dailyCH4(7:16,daily)=CH4V_d(1:10)/24
            
            !  *** .int soil thermal            
            Simu_dailysoilt(1:11,daily)=soilt_d_simu(1:11)/24.
            ! Simu_dailyst(1:11,daily) = testout(1:11)
            Simu_dailyice(1:10,daily)=ice_d_simu(1:10)/24.
                                                        
            Simu_dailywatertable(1,daily)=zwt_d/24.
            Simu_snowdepth(1,daily)=snow_dsim             
            ! endif
            ! if(yr.ge.(yrlim-first_year+1) .and. days.ge.dylim) goto 650
            ! write (122,1202) zwt,snow_dsim
        enddo ! end of idays
                
        storage=accumulation
        stor_use=Storage/times_storage_use
        if(yr.eq.yrs_eq+yr_length .and. do_co2_da.eq.1)then
            write(*,*)yr,LAI,gpp_yr,NPP_yr,pheno
            write(61,601)year,LAI,gpp_yr,NPP_yr,real(pheno)
        endif      
        ! if(MCMC.ne.1) then
        if (do_co2_da.ne.1) then            
            write(*,*)year,LAI,gpp_yr,NPP_yr,pheno,pheno
            write(61,601)year,LAI,gpp_yr,NPP_yr,Ra_yr,Rh_yr, &
            &   ET,rain_yr,transp_yr,evap_yr,runoff_yr,GL_yr,    &
            &   GW_yr,GR_yr,Pool1,Pool2,Pool3,Pool4,Pool5,   &
            &   Pool6,Pool7,Pool8,out1_yr,out2_yr,out3_yr,   &
            &   out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
        endif            
601     format(i7,",",29(f15.4,","))
        accumulation=0.0
        onset=0
    enddo            !end of simulations multiple years
         
    ! if(MCMC.ne.1)then
    if (do_co2_da.ne.1) then
        first_year = year_seq(1)
        i=1
        do nyear=first_year,first_year+yr_length-1
            if(MOD(nyear,4).eq.0)then
                idays=366
            else
                idays=365
            endif
            idays = 365
            do idayOfnyear=1,idays
                write(62,602)i,nyear,idayOfnyear,(Simu_dailyflux(j,i),j=1,12)
                write(662,6602)i,nyear,idayOfnyear,(Simu_dailyflux14(j,i),j=1,14)
                i=i+1
            enddo
        enddo
        do i=1,daily
            ! write(662,6602)i,year,(Simu_dailyflux14(j,i),j=1,14)
            write(63,603)i,(Simu_dailywater(j,i),j=1,31)
            write(64,604)i,(Simu_dailyCH4(j,i),j=1,16)

            write(65,605)i,(Simu_dailysoilt(j,i),j=1,11)
            write(66,606)i,(Simu_dailyice(j,i),j=1,10)
            write(67,607)i,(Simu_dailywatertable(j,i),j=1,1)
            write(68,608)i,(Simu_snowdepth(j,i),j=1,1)
        enddo
       
602      format(3(i7),",",11(f15.4,","),(f15.4))
6602     format(3(i7,","),13(f15.4,","),(f15.4)) 
603      format((i7),",",30(f15.4,","),(f15.4))
604      format((i7),",",15(f15.4,","),(f15.4))
605      format((i7),",",10(f15.4,","),(f15.4))
606      format((i7),",",9(f15.4,","),(f15.4))
607      format((i7),",",(f15.4))
608      format((i7),",",(f15.4))

    endif 
999      continue
    return
        enddo
    end subroutine teco_simu

end module driver