module mod_soil
    use mod_data

    implicit none
    real FLDCAP

    contains
    ! subroutine for soil moisture
    subroutine soilwater()
        real :: infilt_max = 15.
        real DWCL(10), Tr_ratio(10)
        real SRDT(10), depth(10), rain_new, rain_t
        integer nfr
        real infilt_dbmemo, twtadd, wtadd, omegaL(10)
        real exchangeL,supply,demand
        real Tsrdt, tr_allo
        real plantup(10), vtot
        real zmax,thetasmin,zthetasmin,az
        real zwt1,zwt2,zwt3
        real fw(10), ome(10)

        WILTPT = wsmin/100.000
        FLDCAP = wsmax/100.000
        
        do i = 1,10
            dwcl(i)  = 0.0
            evapl(i) = 0.0
            WUPL(i)  = 0.0
            SRDT(i)  = 0.0
            DEPTH(i) = 0.0
        enddo
        ! Layer volume (cm3)
        DEPTH(1) = THKSL(1) ! Determine which layers are reached by the root system. 
        DO i=2,10
            DEPTH(i)=DEPTH(i-1)+THKSL(i)
        enddo
        do i=1,10
            IF(rdepth.GT.DEPTH(i)) nfr=i+1
        enddo
        IF (nfr.GT.10) nfr=10
 
        ! ---------- added for soil thermal    
        if (do_soilphy) then 
            rain_new = rain
            if (ta .lt. -4.) rain_new =0.    ! if (ta .lt. -0.4) rain_new =0.       !dbice   !tuneice     
            rain_t = melt/24+rain_new        ! here it defines how the melt water is added to water input; add melted water hourly
            infilt = infilt+rain_t
            if (ice(1) .gt. 0.0) then        ! Jian: no use in ice? 2022/11/15
                !infilt = 0.0
            endif
        else
            infilt = infilt + rain           ! ..int commented lines for soil thermal module, included in the previous loop; !mm/hour  
        endif     
        
        infilt_dbmemo = infilt    
        ! water infiltration through layers; Loop over all soil layers.
        TWTADD=0
        IF(infilt.GE.0.0)THEN
            ! Add water to this layer, pass extra water to the next.
            WTADD  = AMIN1(INFILT,infilt_max,AMAX1((FLDCAP-wcl(1))*thksl(1)*10.0,0.0)) ! from cm to mm
            ! change water content of this layer
            WCL(1) = (WCL(1)*(thksl(1)*10.0)+WTADD)/(thksl(1)*10.0)
            TWTADD = TWTADD+WTADD       !calculating total added water to soil layers (mm)
            INFILT = INFILT-WTADD !update infilt
        ENDIF
        
        if (do_soilphy) then 
            runoff= INFILT*0.005   !(infilt_rate = 0.0017 defined earlier by Yuan, changed  to 0.001 by shuang )
        else
            runoff=INFILT*0.001    ! Shuang added this elseif line! Shuang Modifed  Mar16 used to be 0.0019, the water lose too much lowest wt was >400
        endif   
        infilt = infilt-runoff
        !----------------------------------------------------------------------------------------------------
        if (transp .gt. 0.2 .and. transp .le. 0.22) then
            infilt = infilt+transp*0.4
        else if (transp .gt. 0.22) then
            ! infilt = infilt+infilt*0.0165
            infilt = infilt+transp*0.8
            ! infilt = infilt+0.22*0.4+(transp-0.22)*0.9
        else
            infilt = infilt+transp*0.001
        endif

        if (evap .ge. 0.1 .and. evap .le. 0.15) then
            infilt = infilt+evap*0.4
        else if (evap .gt. 0.15) then
            infilt = infilt+evap*0.8
        else
            infilt = infilt+evap*0.001
        endif
        !----------------------------------------------------------------------------------------------------   
        ! water redistribution among soil layers
        do i=1,10
            wsc(i) = Amax1(0.00,(wcl(i)-wiltpt)*THKSL(i)*10.0)
            if (do_soilphy) then ! ..int commented lines for soil thermal 
                omegaL(i)=Amax1(0.001,(liq_water(i)*100./thksl(i)-WILTPT)/(FLDCAP-WILTPT))
            else
                omegaL(i)=Amax1(0.001,(wcl(i)-WILTPT)/(FLDCAP-WILTPT))
            endif        
        enddo
        supply = 0.0
        demand = 0.0
        do i=1,9
            if(omegaL(i).gt.0.3)then
                supply    = wsc(i)*(omegaL(i)-0.3)   ! supply=wsc(i)*omegaL(i)
                demand    = (FLDCAP-wcl(i+1))*THKSL(i+1)*10.0*(1.0-omegaL(i+1))
                exchangeL = AMIN1(supply,demand)
                wsc(i)    = wsc(i)- exchangeL
                wsc(i+1)  = wsc(i+1)+ exchangeL
                wcl(i)    = wsc(i)/(THKSL(i)*10.0)+wiltpt
                wcl(i+1)  = wsc(i+1)/(THKSL(i+1)*10.0)+wiltpt
            endif
        enddo
        wsc(10) = wsc(10)-wsc(10)*0.00001     ! Shuang modifed
        runoff  = runoff+wsc(10)*0.00001     ! Shuang modifed
        wcl(10) = wsc(10)/(THKSL(10)*10.0)+wiltpt
        ! end of water redistribution among soil layers
        ! Redistribute evaporation among soil layers
        Tsrdt = 0.0
        DO i=1,10
            ! Fraction of SEVAP supplied by each soil layer
            SRDT(I) = EXP(-6.73*(DEPTH(I)-THKSL(I)/2.0)/100.0) !/1.987 ! SRDT(I)=AMAX1(0.0,SRDT(I)*(wcl(i)-wiltpt)) !*THKSL(I))
            Tsrdt   = Tsrdt+SRDT(i)  ! to normalize SRDT(i)
        enddo
        do i=1,10
            EVAPL(I) = Amax1(AMIN1(evap*SRDT(i)/Tsrdt,wsc(i)),0.0)  !mm
            DWCL(I)  = EVAPL(I)/(THKSL(I)*10.0) !ratio
            wcl(i)   = wcl(i)-DWCL(i)
        enddo
        evap = 0.0       
        do i=1,10
            evap=evap+EVAPL(I)
        enddo
        ! Redistribute transpiration according to root biomass
        ! and available water in each layer
        tr_allo=0.0
        do i=1,nfr
            tr_ratio(i) = FRLEN(i)*wsc(i) !*(wcl(i)-wiltpt)) !*THKSL(I))
            tr_allo     = tr_allo+tr_ratio(i)
        enddo
        do i=1,nfr
            plantup(i) = AMIN1(transp*tr_ratio(i)/tr_allo, wsc(i)) !mm              
            wupl(i)    = plantup(i)/(thksl(i)*10.0)
            wcl(i)     = wcl(i)-wupl(i)
        enddo
        transp = 0.0
        do i=1,nfr
            transp=transp+plantup(i)
        enddo

        ! ---------------------------------------------------------------------------    
        ! water table module starts here
        ! vtot = MAX(145.,wsc(1)+wsc(2)+wsc(3)+infilt)!+wsc(4)+wsc(5)   !total amount of water in top 500mm of soil  mm3/mm2 infilt here is standing water   infilt has distributed to wsc?
        if (do_soilphy) then
            ! vtot = wsc(1)+wsc(2)+wsc(3)+infilt+ice(1)*1000.*(10./9.)+ice(2)*1000.*(10./9.)+ice(3)*1000.*(10./9.)
            ! vtot = wsc(1)+wsc(2)+wsc(3)+ice(1)*1000.*(10./9.)+ice(2)*1000.*(10./9.)+ice(3)*1000.*(10./9.)
            ! vtot = wsc(1)+wsc(2)+wsc(3)+infilt
            vtot = (liq_water(1)+liq_water(2)+liq_water(3))*1000+(ice(1)+ice(2)+ice(3))*1000+infilt
            ! vtot = wsc(1)+wsc(2)+wsc(3)+infilt+ice(1)*1000.*(9./10.)+ice(2)*1000.*(9./10.)+ice(3)*1000.*(9./10.)
            ! write(*,*) ice(1)*1000.,ice(2)*1000.,ice(3)*1000.,wsc(1),liq_water(1)*1000.
        else 
            vtot = wsc(1)+wsc(2)+wsc(3)+infilt
        endif
        ! infilt means standing water according to jiangjiang
        ! vtot = MAX(145.,vtot+145.+rain-evap-transp-runoff)         ! vtot should not be smaller than 145, which is the water content when wt is at -300mm
        phi        = 0.56                           ! soil porosity   mm3/mm3   the same unit with theta
        zmax       = 300                            ! maximum water table depth   mm
        thetasmin  = 0.25                           ! minimum volumetric water content at the soil surface   cm3/cm3
        zthetasmin = 100                            ! maximum depth where evaporation influences soil moisture   mm
        az         = (phi-thetasmin)/zthetasmin     ! gradient in soil moisture resulting from evaporation at the soil surface    mm-1
        
        zwt1 = -sqrt(3.0*(phi*zmax-vtot)/(2.0*az))
        zwt2 = -(3.0*(phi*zmax-vtot)/(2.0*(phi-thetasmin)))
        zwt3 = vtot-phi*zmax                                   
        if ((zwt1 .ge. -100) .and. (zwt1 .le. 0))   zwt = zwt1  !the non-linear part of the water table changing line
        if (zwt2 .lt. -100)                         zwt = zwt2  !the linear part of the water table changing line
        ! if ((zwt2 .lt. -100) .and. (zwt2 .ge. -300))zwt = zwt2 !the linear part of the water table changing line valid when Vtot>145mm
        ! if (zwt2 .le. -300)                         zwt = -300
        if (phi*zmax .lt. vtot)                     zwt = zwt3  !the linear part when the water table is above the soil surface            
        ! water table module ends here
        ! ---------------------------------------------------------------------------------------------------------------
        ! Output fwsoil, omega, and topfws
        ! ..int commented lines below for soil thermal module
        ! do i=1,nfr       
        !     ome(i)=(wcl(i)-WILTPT)/(FLDCAP-WILTPT)
        !     ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
        !     fw(i)=amin1(1.0,3.333*ome(i))
        ! enddo
        ! topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))    
        ! ..int new lines added for soil thermal module 
        do i=1,nfr       
            if (do_soilphy) then 
                ome(i)=(liq_water(i)*100./thksl(i)-WILTPT)/(FLDCAP-WILTPT)
            else 
                ome(i)=(wcl(i)-WILTPT)/(FLDCAP-WILTPT)
                ome(i)=AMIN1(1.0,AMAX1(0.0,ome(i)))
            endif 
            fw(i)=amin1(1.0,3.333*ome(i))
        enddo
        
        if (do_soilphy) then 
            topfws=amax1(0.0,topfws)
        else 
            topfws=amin1(1.0,(wcl(1)-WILTPT)/((FLDCAP-WILTPT)))
        endif     

        fwsoil = 0.0
        omega  = 0.0
        do i=1,nfr
            fwsoil = fwsoil+fw(i)*frlen(i)
            omega  = omega+ome(i)*frlen(i)
        enddo
        return
    end subroutine soilwater

    ! subroutine snow_d(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)
    subroutine snow_d()
        ! real lat,tr,daylength,dec,melt,fa,sublim,dsnow,snow_in,decay_m,fsub
        real tr,daylength,dec,sublim,dsnow, in_snow
        ! real rain_d,snow_dsim,rho_snow,dcount,ta
        ! integer days
        real snow_dsim_pre

        tr        = 0.0174532925
        dec       = sin(((real(iday)-70.)/365.)*360.*tr)*23.44  ! dec=sin(((real(days)-70.)/365.)*360.*tr)*23.44
        daylength = acos(-tan(lat*tr)*tan(dec*tr))/7.5 
        daylength = daylength/tr/24.
                
        if (snow_dsim .ge. 0.) then
            dcount = dcount +1.
        else 
            dcount =0.
        endif
        sublim=0.
        if (ta .gt. 0. .and. snow_dsim .gt. 0.) sublim=fsub*715.5*daylength*esat(ta)/(ta+273.2)*0.001   ! 0.001 from Pa to kPa
        melt=0.
        ! if (ta .gt. 0. .and. snow_dsim .gt. 0.) melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)       !yy version
        if (ta .gt. 1.0e-10 .and. snow_dsim .gt. 0.) melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)   !dbmemo updated version
        ! write(*,*) 'melt=fa*(2.63+2.55*ta+0.0912*ta*rain_d)','fa',fa,'ta',ta,'rain_d',rain_d
        ! if (ta .gt. 0. .and. snow_dsim .gt. 0.) melt=fa*(0.55*ta)     
        if (dcount .gt.0. .and. ta .lt.5.) then
            ! write(*,*)'melt_befor',melt         !dbmemo dbice
            melt=melt*EXP(-decay_m*dcount/365.)  !dbmemo dbice
            ! write(*,*)'melt_after',melt
        endif
        ! write(*,*),EXP(-3.*dcount/365.)
        ! melt=AMIN1(melt, snow_dsim*rho_snow-sublim)
        ! if (melt .lt. 2.) melt=0.
        if (ta .le. 0.) then         ! dbmemo second bug in dbmemo
            in_snow =rain_d
        else
            in_snow = 0.
        endif
        dsnow         = in_snow-sublim-melt 
        snow_dsim_pre = snow_dsim
        snow_dsim     = snow_dsim + dsnow/rho_snow 
        if (snow_dsim .le. 0.0) then 
            snow_dsim=0.0 
            melt = snow_dsim_pre*rho_snow +in_snow-sublim    !! for water part
        endif 
        melt=AMAX1(melt, 0.)
        return
    end subroutine snow_d


    subroutine Tsoil_simu()          
        implicit none 
        real difsv2,difsv1
        real delta
        !real thksl(10)
        real ufw(10),frac_ice1,frac_ice2
        ! real,dimension(10):: Tsoill,liq_water
        real,dimension(11)::testout
        real temph1,temph2     
        real thkns1,thkns2
        real flux_snow
        real condu_water,shcap_water,shcap_ice              !,shcap_snow,condu_snow,depth_ex
        real albedo_water,ice_incr,heat_excess,heat_adjust  !,ice_tw,water_tw
        real inter_var,latent_heat_fusion                   !,QLsoil,Rsoilab3
        real resdh,dnr,dsh,dgh,dle,drsdh                    ! rhocp,slope,psyc,Cmolar,fw1,resoil,rLAI,resht,
        ! real f_om,theta_sat_om,b_om,b_min,phi_om,phi_min,theta_sat,b_tot,phi_sat,gravi
        real water_table_depth,temph_water,temph_snow
        real condu_air,shcap_air,condu(10), shcap(10), condu_ice,tsoill_pre, thd_t
        real ice_density,condu_soil,shcap_soil 
        real resht_lai,snow_depth_t
        real condu_s,tsoill_0,diff_air,d_cor                !,condu_b,dcount,dcount_soil
        real sftmp_pre
        integer n_layers
        real, allocatable ::depth_z(:) 
        n_layers=10
        allocate(depth_z(n_layers))      
        ! write(*,*),thd_t
        ! soil thermal conductivity W m-2 K-1
        ice_density = 916.!916.
        thkns1      = thksl(1)/4.             ! thkns1=thksl(1)/2.
        shcap_ice   = 2117.27*ice_density
        condu_ice   = 2.29
        condu_water = 0.56!0.56
        shcap_water = 4188000.
        condu_soil  = 0.25
        shcap_soil  = 2600000.
        condu_s     = 0.25
        ! thd_t=0.0
        thd_t       = -1.0        
        diff_snow   = 3600.*condu_snow/shcap_snow*10000.
        diff_s      = 3600.*condu_b/shcap_soil*10000.
        latent_heat_fusion = 333700.   ! j kg-1
        condu_air    = 0.023
        shcap_air    = 1255.8
        diff_air     = 3600.*condu_air/shcap_air*10000.      
        water_tw     = zwt*0.001-ice_tw ! might means total water that is liquid, add up all layers
        water_table_depth=zwt*0.1
        snow_depth_t = snow_depth - 0.46*0.0     ! warming in Tair impact on snow_depth
                                                   ! in unit cm 0.46 based on snow_depth vs. tair regression     
        if (snow_depth_t .lt. thd_snow_depth) snow_depth_t =0.0
        if (snow_depth_t .gt. 0.) then
            dcount_soil = dcount_soil +1./24.
        else 
            dcount_soil =0.
        endif
    
        if (water_table_depth .lt. 4. .and. water_table_depth .gt. 0.0) water_table_depth =0.    ! avoid numerical issues when 
        ! if (water_table_depth .lt. -99.) water_table_depth =-30.    ! temporary for NaN    
        ! if (water_table_depth .lt. -299.) water_table_depth =-30.    ! -299 is a more reasonable value Shuang Ma         
        albedo_water = 0.1      
        ! soil water conditions
        WILTPT       = wsmin/100.
        FILDCP       = wsmax/100.
        TairK        = Tair+273.2            
        flux_snow    = 0.0   
        depth_z      = (/0., 0., 0., 0., 0., 0., 0.,0.,0.,0./) 
        ! ..int add unfrozen water ratio
        ! ufw=(/0.0042,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063,0.0063/)
        ufw=(/0.0163,0.0263,0.0563,0.0563,0.0563,0.1162,0.1162,0.1162,0.1162,0.1162/)
        ! ufw=(/0.0042,0.009,0.009,0.0563,0.0563,0.1162,0.1162,0.1162,0.1162,0.1162/)
        frac_ice1 = 0.01!0.015
        frac_ice2 = 0.001!0.01
        ! if (snow_depth_t .gt. 0.0) then 
        !     emsoil =0.98
        ! elseif (water_table_depth .gt. 0.0) then
        !     emsoil =0.99
        ! endif
        QLsoil   = emsoil*sigma*((sftmp+273.2)**4)
        Rsoilab3 = (QLair+QLleaf)*(1.0-rhoS(3))-QLsoil       
        ! Total radiation absorbed by soil
        if (snow_depth_t .gt. 0.0) then 
            Rsoilabs = (Rsoilab1+Rsoilab2)*(1-albedo_snow)/(1-0.1)+Rsoilab3  
        elseif (water_table_depth .gt. 0.0) then 
            Rsoilabs = (Rsoilab1+Rsoilab2)*(1-albedo_water)/(1-0.1)+Rsoilab3  
        else
            Rsoilabs = Rsoilab1+Rsoilab2+Rsoilab3
        endif
          
        ! thermodynamic parameters for air
        rhocp  = cpair*Patm*AirMa/(Rconst*TairK)      
        H2OLv  = H2oLv0-2.365e3*Tair
        slope  = (esat(Tair+0.01)-esat(Tair))/0.01   
    
        psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar = Patm/(Rconst*TairK)
        fw1    = AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.3),1.0)    

        if (water_table_depth .gt. 0.0) then 
            Rsoil = 0. 
        else 
            Rsoil=30.*exp(0.2/fw1)
        endif 
        ! Rsoil=40.
        ! Rsoil=5.
        rLAI=exp(FLAIT)     
        ! latent heat flux into air from soil
        !       Eleaf(ileaf)=1.0*
        ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
        ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
    
        Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
               &      (slope+psyc*(Rsoil/(raero+rLAI)+1.))
        resht_lai=resht*FLAIT
        ! resht_lai= resht*exp(FLAIT)/15. ! need improvement, should be a function of LAI 
        ! if (water_table_depth .gt. 0.0) resht_lai=resht/FLAIT*0.2

        ! resht_lai=200.      
        Hsoil=rhocp*(sftmp-Tair)/resht_lai   
        ! Hsoil=1010.*1.17*(sftmp-Tair)/resht_lai
        i=1;
        condu(i) = (FILDCP-wcl(i))*condu_air+liq_water(i)/(thksl(i)*0.01)*condu_water+ &
                    &  ice(i)/(thksl(i)*0.01)*condu_ice +(1-FILDCP)*condu_soil
        shcap(i) = (FILDCP-wcl(i))*shcap_air+liq_water(i)/(thksl(i)*0.01)*shcap_water+ &
                    &  ice(i)/(thksl(i)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil
        difsv1   = 3600.*condu(i)/shcap(i)*10000.    
        G        = condu(1)*(sftmp-tsoill(1))/(thksl(1)/2.*0.01)
        if (snow_depth_t .gt. 0.0) then 
            G = condu_snow*(sftmp-Tsnow)/(snow_depth_t/2.*0.01)
        endif
          
        ! thksl(1)
        ! G=0.   
        ! Residual heat energy.
        RESDH = Rsoilabs-Hsoil-Esoil-G
        ! First derivative of net radiation; sensible heat; ground heat;
        DNR   = 4.*emsoil*sigma*(sftmp+273.2)**3
        DSH   = rhocp/resht_lai 
        DGH   = condu_s/(thksl(1)/2.*0.01)
        DLE   = (DNR+DGH)*slope/(slope+psyc*(Rsoil/(raero+rLAI)+1.))      
        drsdh = -DNR-DSH-DGH-DLE
        ! Calculate increment DELTA.
        DELTA     = resdh/drsdh
        sftmp_pre = sftmp
        sftmp     = sftmp-DELTA
        if (ABS(sftmp_pre -sftmp) .gt. 20. ) sftmp=sftmp_pre  
        tsoill_0  = sftmp

        do i=1,10
            Tsoill_pre = tsoill(i)    
            if (water_table_depth .lt. 0.0 .and. -water_table_depth .lt. depth_z(i)) then
                liq_water(i) = FILDCP*thksl(i)*0.01-ice(i)
            else
                liq_water(i) = wcl(i)*thksl(i)*0.01-ice(i)
            endif          
            if (i .eq. 1) then 
                depth_z(1) = thksl(1)
            else 
                depth_z(i) = depth_z(i-1)+thksl(i)
            endif
            
            thkns2 = (thksl(i)+thksl(i+1))/2.
                     
            if (i .eq. 10) then
                difsv2 = 3600.*condu(i)/shcap(i)*10000. 
            else
                condu(i+1) = (FILDCP-wcl(i+1))*condu_air+liq_water(i+1)/(thksl(i+1)*0.01)*condu_water+ &
                                &  ice(i+1)/(thksl(i+1)*0.01)*condu_ice +(1-FILDCP)*condu_soil
                shcap(i+1) = (FILDCP-wcl(i+1))*shcap_air+liq_water(i+1)/(thksl(i+1)*0.01)*shcap_water+ &
                                &  ice(i+1)/(thksl(i+1)*0.01)*shcap_ice +(1-FILDCP)*shcap_soil 
                difsv2     = 3600.*condu(i+1)/shcap(i+1)*10000.
            endif     
            temph2 = (difsv1+difsv2)*(Tsoill(i)-Tsoill(i+1))/thkns2 
            !!!!!!!!!!!!!!!!!!!! start first layer !!!!!!!!!!!!!!!!!!!!!!
            !!!!!! adjust if there are snow or water layer above !!!!!!!!!!!!!!!!!!!!
            if(i.eq.1) then
                if (snow_depth_t .gt. 0.) then   
                    temph_snow = Amin1(diff_snow,difsv1)*(Tsnow-Tsoill(1))/((snow_depth_t+thksl(1))/2.)
                    Tsnow      = Tsnow+(exp(-depth_ex*snow_depth_t)*diff_snow*(sftmp-Tsnow)/(snow_depth_t/2.) &
                                    &    -temph_snow)/(snow_depth_t/2.+(snow_depth_t+thksl(1))/2.) 
                    Tsoill(1)  = Tsoill(1)+(temph_snow &
                                    &    -temph2)/((snow_depth_t+thksl(1))/2.+thkns2) 
                    if (Tsnow .gt.0.0) then 
                        Tsnow     = 0.0   
                        Tsoill(1) = 0.
                    endif
                    drsdh    = 0.0    ! temporarily set drsdh =0 for heat adjustment of soil when  
                    tsoill_0 = (Tsoill(1)+Tsnow)/2.
                elseif (water_table_depth .gt. 0.) then  
                    temph_water = (3600.*condu_water/shcap_water*10000.+difsv1)*(Twater-Tsoill(1))/((water_table_depth+thksl(1))/2.)! there is snow layer 
                    Twater      = Twater+(2.*3600.*condu_water/shcap_water*10000.*(sftmp-Twater)/(water_table_depth/2.) &
                                    &        -temph_water)/(water_table_depth/2.+(water_table_depth+thksl(1))/2.) 
                    !!!!!!!!!!!!!!!!!!  Phase change surface water !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if (Twater .lt. 0.0 .and. water_tw .gt. 0.0) then  ! freeze 
                        heat_excess = -(shcap_water/360000.*water_tw*100.-drsdh)*Twater
                        ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density
                        !!-----------int add mechanism of unfrozen water in frozen soil layers, typically happens in high latitude region
                        !!              according to obs soil water content, winter water never goes below 0.063 at -20cm and 0.042 at surface layer
                        !!       !tuneice
                        !                if (ice_incr .lt. 0.) then   
                        !                   if (i .eq. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = 0. !ice_incr*0.1
                        !                          ice_incr = ice_incr*frac_ice1
                        !                       endif
                        !    !               elseif (i .eq. 2) then
                        !    !                   if (liq_water(i) .le. 0.063) then
                        !    !                       ice_incr = 0.
                        !    !                   endif
                        !                   elseif ( i .gt. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = ice_incr*0. !0.9
                        !                           ice_incr = ice_incr*frac_ice2
                        !                       endif
                        !                   endif
                        !                endif
                        !! $$$$$$$$$$$   
                       
                        ! write(*,*)'water_tw',water_tw
                        if (ice_incr .lt. water_tw) then
                            ice_tw   = ice_tw +ice_incr
                            water_tw = water_tw-ice_incr
                            Twater   = 0.0
                            Tice     = 0.0
                        else
                            ice_tw   = ice_tw +water_tw
                            water_tw = 0.0
                            Tice     = Tice - latent_heat_fusion*(ice_incr-water_tw)*ice_density/(shcap_ice*ice_tw)
                        endif     
                    elseif (Twater .gt. 0.0 .and. ice_tw .gt. 0.0) then    ! thraw              
                        heat_excess  = (shcap_water/360000.*ice_tw*100.-drsdh)*Twater
                        ice_incr     = heat_excess*3600./latent_heat_fusion/ice_density
                        !! $$$$$$$$$$$
                        !! $$$$$$$$$$$      !tuneice
                        !                if (ice_incr .lt. 0.) then   
                        !                   if (i .eq. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = 0. !ice_incr*0.1
                        !                          ice_incr = ice_incr*frac_ice1
                        !                       endif
                        !    !               elseif (i .eq. 2) then
                        !    !                   if (liq_water(i) .le. 0.063) then
                        !    !                       ice_incr = 0.
                        !    !                   endif
                        !                   elseif ( i .gt. 1.) then
                        !                       if (liq_water(i) .le. ufw(i)) then
                        !!                           ice_incr = ice_incr*0. !0.9
                        !                           ice_incr = ice_incr*frac_ice2
                        !                       endif
                        !                   endif
                        !                endif
                        !! $$$$$$$$$$$                   
                        !! $$$$$$$$$$$                   
             
                        if (ice_incr .lt. ice_tw) then
                            ice_tw   = ice_tw -ice_incr
                            water_tw = water_tw+ice_incr
                            Twater   = 0.0
                            Tice     = 0.0
                        else
                            water_tw = water_tw +ice_tw
                            ice_tw   = 0.0
                            Twater   = Twater + latent_heat_fusion*(ice_incr-ice_tw)*ice_density/(shcap_water*water_tw)
                        endif
                    endif                       
                    !!!!!!!!!!!!!!!!!!!!!!!!! end of phase change for surface layer !!!!!!!!!!!!!!!!!!!  
                    temph2=(difsv1+3600.*condu_water/shcap_water*10000.)*(Tsoill(i)-Tsoill(i+1))/thkns2 
                    if (water_tw .eq. 0.0 .and. ice_tw .gt. 0.0) then 
                        Tsoill(1) = Tsoill(1)+(2.*3600.*condu_ice/shcap_ice*10000.*(Tice-Tsoill(1))/thkns1 &
                                    &     -temph2)/(thkns1+thkns2) 
                    else 
                        Tsoill(1) = Tsoill(1)+(2.*3600.*condu_water/shcap_water*10000.*(Twater-Tsoill(1))/thkns1 &
                                    &     -temph2)/(thkns1+thkns2) 
                    endif
                    drsdh = 0.0    ! temporarily set drsdh =0 for heat adjustment of soil       
                else   
                    Tsoill(1) = Tsoill(1)+(diff_s*(sftmp-Tsoill(1))/thkns1 &
                                &     -temph2)/(thkns1+thkns2)
                endif
                !!!!!  phase change in top soil       
                heat_excess = drsdh*(thd_t-Tsoill(i))+shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.         
                ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density        
                !
                !! $$$$$$$$$$$
                !! $$$$$$$$$$$      !tuneice
                !                if (ice_incr .lt. 0.) then   
                !                   if (i .eq. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = 0. !ice_incr*0.1
                !                          ice_incr = ice_incr*frac_ice1
                !                       endif
                !    !               elseif (i .eq. 2) then
                !    !                   if (liq_water(i) .le. 0.063) then
                !    !                       ice_incr = 0.
                !    !                   endif
                !                   elseif ( i .gt. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = ice_incr*0. !0.9
                !                           ice_incr = ice_incr*frac_ice2
                !                       endif
                !                   endif
                !                endif
                !! $$$$$$$$$$$   
                !! $$$$$$$$$$$                           
                !!          
                inter_var = ice(i)   
                if (ice_incr .lt. 0.) then     ! freeze             
                    ice(i) = Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)            
                else 
                    ice(i) = Amax1(ice(i)-ice_incr,0.0)              
                endif
                !! readjust energy and temp 
                heat_adjust = heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
                Tsoill(i)   = thd_t+heat_adjust/(shcap(i)*thksl(i)/360000.-drsdh)      
            else
                ! if ( i .gt. 9) then 
                !     temph2=0
                !     thkns2=500  ! boundary conditions, rethink
                ! endif
                if ( i .gt. 9) then 
                    temph2 = 0.00003
                    thkns2 = 500  ! boundary conditions, rethink
                endif            
                Tsoill(i)   = Tsoill(i)+(temph1-temph2)/(thkns1+thkns2)    
                heat_excess = shcap(i)*thksl(i)*(Tsoill(i)-thd_t)/360000.        
                ice_incr    = heat_excess*3600./latent_heat_fusion/ice_density         
                !! $$$$$$$$$$$
                !! $$$$$$$$$$$      !tuneice
                !                if (ice_incr .lt. 0.) then   
                !                   if (i .eq. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = 0. !ice_incr*0.1
                !                          ice_incr = ice_incr*frac_ice1
                !                       endif
                !    !               elseif (i .eq. 2) then
                !    !                   if (liq_water(i) .le. 0.063) then
                !    !                       ice_incr = 0.
                !    !                   endif
                !                   elseif ( i .gt. 1.) then
                !                       if (liq_water(i) .le. ufw(i)) then
                !!                           ice_incr = ice_incr*0. !0.9
                !                           ice_incr = ice_incr*frac_ice2
                !                       endif
                !                   endif
                !                endif
                !! $$$$$$$$$$$   
                !! $$$$$$$$$$$                   
                !            
                
                inter_var = ice(i) 
                if (ice_incr .lt. 0.) then     ! freeze             
                    ice(i) = Amin1(liq_water(i)+inter_var,ice(i)-ice_incr)             
                else 
                    ice(i) = Amax1(ice(i)-ice_incr,0.0)              
                endif         
                !! readjust energy and temp 
                heat_adjust = heat_excess-latent_heat_fusion*(inter_var-ice(i))*ice_density/3600.
                Tsoill(i)   = thd_t+heat_adjust/(shcap(i)/360000.*thksl(i))
            endif
           
            if (ABS(tsoill_pre -tsoill(i)) .gt. 5. ) Tsoill(i)=tsoill_pre
            TEMPH1 = TEMPH2
            THKNS1 = THKNS2
            DIFSV1 = DIFSV2        
        enddo
        testout(1)    = tsoill_0       ! testout(1)=tsnow
        testout(2:11) = tsoill(1:10) 
        deallocate(depth_z)
        return 
    end subroutine Tsoil_simu
    
    subroutine methane(Rh_pools,Tsoil,zwt,wsc,      &      !update single value in a hourly loop when MEMCMC=0
        &             phi,LAIMIN,LAIMAX,           &
        &             ProCH4,Pro_sum,OxiCH4,Oxi_sum,Fdifu,Ebu_sum,Pla_sum,simuCH4,CH4,CH4_V,   &
        ! &             ProCH4,Pro_sum,OxiCH4,Oxi_sum,Fdifu,Ebu_sum,Pla_sum,simuCH4,CH4,  &
        &             r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi, &
        &             testout, do_soilphy)  !update single value of Rh_pools,Tsoil,zwt,wsc 
        !                                                           in a hourly loop when MEMCMC=1
        !******************************************************************************************************************************
        !****!introduce variables and constants used in this subroutine
        !****************************************************************************************************************************** 
        !     set soil layers
        !****************************************************************************************************************************** 
        implicit none
        ! integer i,MEMCMC
        integer i
        integer,parameter :: nlayers=10       !use this statement to set the parameter value
        real zwt    
        real consum
        !****************************************************************************************************************************** 
        !     set values for MEMCMC
        !******************************************************************************************************************************       
        integer,parameter :: miterms=17
        integer,parameter :: ilines=9000      
        !******************************************************************************************************************************
        !     CH4 Production      
        !******************************************************************************************************************************      
        !      real Rhetero
        real Rh(nlayers),Rh_pools(5),Rh_h,ProCH4(nlayers),Pro_sum
        real r_me         !release ratio of CH4 to CO2
        real Q10pro
        real fSTP(nlayers)         !CH4 production factor of soil temperature
        real vt,xt
        real Tmax_me,Tpro_me
        real fpH          !CH4 production factor of soil pH
        real fEhP         !CH4 production factor of soil redox potential
        real FRLEN(nlayers)        !fraction of root in each layer
        real Tsoil
        !****************************************************************************************************************************** 
        !     CH4 Oxidation
        !******************************************************************************************************************************      
        !      real CH4(nlayers),CH4_V(nlayers+1)          !both are CH4 concentration: CH4(nlayers)unit gC/m2, CH4_V(nlayers) unit g C/ m3
        real CH4(nlayers),CH4_V(nlayers)          !both are CH4 concentration: CH4(nlayers)unit gC/m2, CH4_V(nlayers) unit g C/ m3
        real wsc(nlayers)      
        real OxiCH4(nlayers),Oxi_sum       !CH4 oxidation
        real Omax_layers(nlayers),Omax       !maximum oxidation rate
        real kCH4_layers(nlayers),kCH4       !system specific half saturation constant
        real Q10oxi
        real fCH4(nlayers)         !CH4 oxidation factor of CH4 concentration
        real fSTO(nlayers)         !CH4 oxidation factor of soil temperature
        real fEhO         !CH4 oxidation factor of soil redox potential
        real Toxi
        !****************************************************************************************************************************** 
        !     CH4 Diffusion
        !******************************************************************************************************************************      
        real Deff(nlayers)     !CH4 effective diffusivity !!used in mineral soil  v1.1 
        real D_CH4_a           !CH4 diffusion coefficient in air  !unit cm2 s-1   diffusivity of CH4 in air
        real D_CH4_w           !CH4 diffusion coefficient in water  !unit cm2 s-1   diffusivity of CH4 in water
        real phi          !soil porosity  also used in water table part
        real fwater(nlayers),fair(nlayers)
        real D_CH4_soil(nlayers),D_CH4_soil_a(nlayers),D_CH4_soil_b(nlayers)      !!used in organic peat soil  v1.2
        real fcoarse      !relative volume of coarse pores depending on soil texture  Zhuang 2004
        real ftort        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
        !suggesting that the distance covered by diffusion is about two thirds of the length of the real average path
        real SAND         !relative contents of sand (%) in the soil
        real PVSAND       !relative volume of coarse pores in sandy soils     set to 0.45     value from Walter 2001
        real SILT         !relative contents of silt (%) in the soil
        real PVSILT       !relative volume of coarse pores in silty soils     set to 0.20     value from Walter 2001
        real CLAY         !relative contents of clay (%) in the soil
        real PVCLAY       !relative volume of coarse pores in clayish soils     set to 0.14   value from Walter 2001
        real DEPTH(10)        !depth in soil  will define it inside this subroutine again      resolution 100mm 200mm
        real THKSL(10)        !will define it inside this subroutine again  
        ! real Fdifu(nlayers+1)
        real Fdifu(nlayers)
        real CH4_atm      !concentration of CH4 in atmosphere     seen as 0 cause the value is too low someone use 0.076
        real simuCH4      !simulated CH4 emission
        ! ***********  Boundary condition parameters    *************      
        real ScCH4                                 !Schmidt numbers for methane Wania
        real pistonv                               !Piston velocity
        real Ceq                                   !equilibrium concentration of gas in the atmosphere
        real kHinv                                 !Henry's coefficient dependent variable on left side of equation, T is the independent variable
        real kH_CH4         !Henry's constant at standard temperature (CH4) Unit L atm mol-1
        real CHinv          !Coefficient in Henry's Law Unit K      
        real Tsta           !standard temperature Unit K
        real Ppartial       !CH4 partial pressure in air Unit atm
    
        !****************************************************************************************************************************** 
        !     Ebullition 
        !******************************************************************************************************************************      
        real CH4_thre,CH4_thre_ly(nlayers),EbuCH4(nlayers),Kebu
        real Ebu_sum_unsat,Ebu_sum_sat,Ebu_sum          !sum value one dimension is enough 
        integer wtlevelindex
        !****************************************************************************************************************************** 
        !     Plant transport
        !******************************************************************************************************************************      
        real PlaCH4(nlayers),Pla_sum
        real LAIMIN,LAIMAX
        real Tveg,Tgr,Tmat,fgrow,Pox,Kpla
        !****************************************************************************************************************************** 
        !******************************************************************************************************************************
        ! Yuan added for soil temp  
        logical do_soilphy
        real testout(11), tsoil_layer(11)
        !****************************************************************************************************************************** 
        !******************************************************************************************************************************
        !      MEMCMC=0   ! note here, any changes here result unexpected bug 
        Rh_h=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
        tsoil_layer = testout
        FRLEN = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)             
        ! FRLEN = (/0.1,0.25,0.25,0.2,0.1,0.05,0.025,0.015,0.005,0.005/)
        ! FRLEN = (/0.05,0.1,0.1,0.1,0.15,0.25,0.25,0.0,0.00,0.00/)
        thksl = (/10.,10.,10.,10.,10.,20.,20.,20.,20.,20./)
        simuCH4 = 0.0                 ! v1.2 
        do i = 1, nlayers
            !!!!!!!put it out of the subroutine
            !****************************************************
            !* Rh weighed according to the distribution of root *
            !****************************************************
            if (i .LE. 3) then                                 ! the empirical method used here is from CLM4.5
                Rh(i)= 0.5*Rh_h*FRLEN(i)+((0.5*Rh_h)/0.3)*0.1   
                ! Rh(h,i)Rh produced by each layer per hour  unit should be g C m-2 h-1 
            else                                               ! i*10: depth of ith soil layers
                Rh(i)= 0.5*Rh_h*FRLEN(i)
            endif
            ! Rh(i) = Rh(i) + OxiCH4(i)*(11/4)
        enddo   
           
        !****************************************************************************************************************************** 
        !******************************************************************************************************************************            
    
        !****************************************************          
        !A. methane production     hourly  gC m-2 hour-1
        !Methane production is modeled as an anaerobic process that occurs in the saturated zone of the soil profile ZHUANG
        !****************************************************
        !Rh_h=Rh_pools(1)+Rh_pools(2)+Rh_pools(3)+Rh_pools(4)+Rh_pools(5)  !hourly Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
        !r assignment
        ! r_me=0.3      !find in parafile
        Tmax_me=45.0
        ! Tpro_me=10.0
        ! Q10pro=3.0    !find in parafile
        do i = 1,nlayers          
            if (do_soilphy) then
                if (tsoil_layer(i+1) .lt. 0.0) then
                    fSTP(i) = 0.0
                else if (tsoil_layer(i+1) .gt. Tmax_me) then
                    fSTP(i) = 0.0
                else if (tsoil_layer(i+1) .ge. 0.0 .and. tsoil_layer(i) .le. Tmax_me) then
                    fSTP(i) = Q10pro**((tsoil_layer(i+1)-Tpro_me)/10)        !Tsoil is the only variable
                endif
            else 
                if (Tsoil .lt. 0.0) then
                    fSTP(i) = 0.0
                else if (Tsoil .gt. Tmax_me) then
                    fSTP(i) = 0.0
                else if (Tsoil .ge. 0.0 .and. Tsoil .le. Tmax_me) then
                    fSTP(i) = Q10pro**((Tsoil-Tpro_me)/10)        !Tsoil is the only variable
                endif
            endif
        enddo
        ! fpH assignment
        fpH=1.0
        ! fEhP assignment
        fEhP=1.0
        
        depth(1)=10.0                                  !calculate soil depth unit cm
        do i=2,nlayers
            depth(i)=depth(i-1)+THKSL(i)
        enddo
          
        Pro_sum=0.0
        do i = 1,nlayers
            ! (depth(i)*10)                   !convert unit from cm to mm
            ! (THKSL(i)*10)                   !convert unit from cm to mm convert the unit in each of the equations
            if ((depth(i)*10) .le. -zwt) then
                ProCH4(i)=0.0
            else
                if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                    ProCH4(i)=Rh(i)*r_me*fSTP(i)*fpH*fEhP*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))     ! *percent
                elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                    ProCH4(i)=Rh(i)*r_me*fSTP(i)*fpH*fEhP
                endif
            endif
            Pro_sum=Pro_sum+ProCH4(i)
        enddo
    
        !**************************************************
        !Add CH4 production to CH4 pool    (gC layer -1)=(gC m-2)
        !**************************************************
    
        do i=1,nlayers
            CH4(i) = CH4(i) + ProCH4(i)
        enddo  
        ! END OF METHANE PRODUCTION
        ! ********************************************************************************************************************
        ! B. methane oxidation      hourly  unit gC m-2 h-1     !!!!!!!!!!!method of CLM and Zhuang!!!!!!!!!!!!
        ! Methane oxidation is modeled as an aerobic process that occurs in the unsaturated zone of the soil profile ZHUANG
        ! ********************************************************************************************************************
        ! fSTO assignment
        ! ***************          
        Q10oxi=2.0      !Zhu 2014 results from previous studies  unit 1  also used by zhang
        ! Toxi=10.0       !Zhuang 2004 table1 Boreal Forest Wetland
        do i=1,nlayers
            if (do_soilphy) then
                fSTO(i)=Q10oxi**((tsoil_layer(i+1)-Toxi)/10.0)
            else
                fSTO(i)=Q10oxi**((Tsoil-Toxi)/10.0)
            endif
        enddo
        ! fEhO assignment
        fEhO=1.0        !Walter 2000  did not consider it, equal to value of 1
    
        ! Omax assignment
        ! ***************
        Oxi_sum=0.0
        do i = 1,nlayers
            ! Omax=1.5
            ! Omax=15.0  !!find in parafile Zhuang 2004 table1 Boreal Forest Wetland μmol L-1 h-1 system specific maximum oxidation coefficient
            ! convert the unit of Omax from μmol L-1 h-1 to gC m-2 h-1
            ! /1000,000 to get mol
            ! *12 cmass to get gC
            ! *1000 to get from dm-3(L) to m-3
            ! *(wsc*0.001) to get unit of omax_layers from m-3 to m-2     !caution that wsc unit is mm
            ! ** w  /   (w/t)           CLM used 
            Omax_layers(i)=(Omax/(1000000))*12*1000*(wsc(i)*0.001)     !convert the unit of Omax from μmol L-1 h-1 to gC m-2 h-1
            ! Omax_layers(i)=(Omax/(1000000))*12*1000*(THKSL(i)*10.0)*0.001     !modified on 11/27/2016 no sig change in oxidation and emission
            ! in unsaturated part of oxidation in CLM, they used the Omax/10 but did not expained why   they also /water volume
            ! fCH4 assignment
            ! ***************          
            ! kCH4=5.0     !!find in parafile Zhuang 2004 range between 1 and 66.2μmol L-1 system specific half saturation constant  1.0e
            ! convert the unit of kCH4 from μmol L-1 to gC m-2
            kCH4_layers(i)=(kCH4/(1000000))*12*1000*(wsc(i)*0.001)    !convert the unit of kCH4 from μmol L-1 to gC m-2
            ! then calculate fCH4 with CH4(i) and kCH4_layers(i) 
            fCH4(i)=CH4(i)/(kCH4_layers(i)+CH4(i))   !  CH4 concentration factor
            if ((depth(i)*10.0) .le. -zwt) then                !unit of Omax: gC m-2 h-1
                OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO!*0.1      !wrong:*(THKSL(i)/1000)!mm to m account for the thickness
                ! OxiCH4(i)=CH4(i)*0.001
            else
                if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then
                    if (i .eq. 1) then
                        OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO*((-zwt)/(THKSL(i)*10.0))
                    else
                        OxiCH4(i)=Omax_layers(i)*fCH4(i)*fSTO(i)*fEhO*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))      !  *percent
                    endif
                    ! OxiCH4(i)=CH4(i)*0.001*(((-zwt)-(depth(i-1)*10.0))/(THKSL(i)*10.0))
                else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then
                    OxiCH4(i)= 0.0
                endif
            endif            
            if (OxiCH4(i) .gt. CH4(i)) then
                OxiCH4(i)=CH4(i)
            endif  
            Oxi_sum=Oxi_sum+OxiCH4(i)   
        enddo 
    
        !*******************************************************************
        !minus CH4 oxidation from CH4 pool     
        !*******************************************************************
        do i=1,nlayers
            CH4(i) = CH4(i) - OxiCH4(i)               !minus CH4 oxidation from CH4 pool
            CH4_V(i) = CH4(i)/(wsc(i)*0.001)          !convert concentration from gC/m2 to gC/m3
                                                    !CH4_V(i) can be used for DA with observation data in soil layers
        enddo
    
          
        ! END OF METHANE OXIDATION
          
        ! ****************************************************
        ! C. methane diffusion
        ! ****************************************************
        ! Parameters assignment 
        D_CH4_a=0.2            !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in air
        D_CH4_a=(D_CH4_a/10000.0)*3600.0        !unit m2 h-1
        D_CH4_w=0.00002        !unit cm2 s-1   D_CH4_a is the molecular diffusion coefficient of methane in water
        D_CH4_w=(D_CH4_w/10000.0)*3600.0        !unit m2 h-1          
        ftort=0.66        !tortuosity coefficient with a value of 0.66    Walter and Heimann 2000
        ! parameters for fcoarse algorithm      
        SAND=0.4             !   %   SPRUCE site value    0.4
        SILT=0.4             !   %   SPRUCE site value   0.4
        CLAY=0.2             !   %   SPRUCE site value   0.2
        PVSAND=0.45       !relative volume of coarse pores in sandy soils       set to 0.45     value from Walter 2001 zhuang
        PVSILT=0.20       !relative volume of coarse pores in silty soils       set to 0.20     value from Walter 2001 zhuang
        PVCLAY=0.14       !relative volume of coarse pores in clayish soils     set to 0.14     value from Walter 2001 zhuang  
        fcoarse=SAND*PVSAND+SILT*PVSILT+CLAY*PVCLAY
        CH4_atm=0.076       !unit umol L-1
        ! CH4_atm=0.0       !unit umol L-1      
        ! ******************************************************************************************************
        ! * Peat soil solution for diffusion coefficient: Equations for D_CH4_soil *         v1.2    Millington and Quirk Model
        ! ******************************************************************************************************
        do i=1,nlayers
            fwater(i) = wsc(i)/(THKSL(i)*10)      
            fair(i) = phi-fwater(i)
                        
            D_CH4_soil_a(i) = (((fair(i))**(10/3))/((phi)**2))*D_CH4_a
            D_CH4_soil_b(i) = D_CH4_W
            if (fair(i) .ge. 0.05) then
                D_CH4_soil(i) = D_CH4_soil_a(i)
            else
                D_CH4_soil(i) = D_CH4_soil_b(i)
            endif
            ! D_CH4_soil(i) = ge(fair,0.05)*D_CH4_soil_a(i) + lt(fair,0.05)*D_CH4_soil_b(i)        
                    
            ! Here I divided into saturated layer and unsaturated layer conditions because in most cases fair is > 0.05 and that there might be too much diffusion v1.2
            ! or maybe I can adjust the value of threshold 0.05 to around 0.08 as in most cases fwater=0.88 fair=0.07
            ! if (zwt .ge. 0.0) then                                  !when water table is above the soil surface
            !     Deff(i) = D_CH4_W
            ! elseif (zwt .lt. 0.0) then                                  !when water table is below the soil surface
            !     if ((depth(i)*10.0) .le. -zwt) then               !acrotelm layers
            !         Deff(i) = D_CH4_soil(i)
            !     elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then       !partly acrotelm layer
            !         Deff(i) = D_CH4_soil(i)
            !     elseif (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then   !catotelm layers
            !         Deff(i) = D_CH4_W
            !     endif
            ! endif
            ! in this case diffusion should be more
            Deff(i) = D_CH4_soil(i)
        enddo 
    
            
        ! ******************************************************************************************************
        ! * Mineral soil solution for diffusion coefficient: Equations for D_CH4_soil *         v1.1   Three-porosity-model
        ! ******************************************************************************************************
        ! do i = 1,nlayers
        !     fwater = wsc(i)/(THKSL(i)*10)
        !     fair = phi-fwater
        !     fwater = 0.68         ! switch on when testing the effect of fwater on diffusion 0.6 crash 0.7fine  02172017
        !     Deff(i) = D_CH4_a*fcoarse*ftort*phi*(phi-fwater)+D_CH4_w*fwater           !
        ! enddo
    
          
        !convert the unit of CH4_atm from μmol L-1 to gC m-3
        !/1000,000 to get mol
        !*12 cmass to get gC
        !*1000 to get from dm-3(L) to m-3
        CH4_atm = (CH4_atm/1000000)*12*1000
              
        ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! New improvement 2017: Boundary condition     
        kH_CH4 = 714.29
        CHinv = 1600.0
        Tsta = 298.15
        ! Ppartial = 1.7E-6
        Ppartial = 1.7E-20 
          
        ScCH4 = 1898 - 110.1*Tsoil + 2.834*Tsoil**2 - 0.02791*Tsoil**3
        pistonv = 2.07 * (ScCH4/600)**(-1/2)
        kHinv = kH_CH4 /((exp(CHinv*(1/(Tsoil+273.15)-1/Tsta))))
        ! write (*,*) kHinv,pistonv
        Ceq = Ppartial / kHinv    ! Ceq: mol L-1   p_partial: atm  kHinv：L atm mol-1
        ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        Fdifu(1) =  pistonv * (CH4_V(1) - Ceq)
        ! Fdifu(1) =  -pistonv * (CH4_V(1) - Ceq)
        ! if (zwt .ge. -100.0) then
        ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! Fdifu(1) = - pistonv * (CH4_V(1) - Ceq)                         !switch on/off
        ! else
        ! Fdifu(1) = Deff(1)*(CH4_V(1)-CH4_atm)/(THKSL(1)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! endif
             
        ! if (zwt .ge. -200 .and. zwt .le. 100.0) then
        ! Fdifu(2) = - pistonv * (CH4_V(2) - Ceq) 
        ! else
        ! Fdifu(2) = Deff(2)*(CH4_V(2)-CH4_V(1))/(THKSL(2)*0.01)         !refer to the interface of methane flux from layer 1 to atmosphere  cm to m   switch on/off
        ! endif   
    
        do i = 2,nlayers                                  !refer to flux from layer ii to ii-1 
            Fdifu(i)= Deff(i)*(CH4_V(i)-CH4_V(i-1))/(THKSL(i)*0.01)      !the unit of Fdifu is gC/m-2/h
        enddo
        ! CH4_V(11) = CH4_V(10)
        ! Fdifu(11) = Deff(10)*(CH4_V(11)-CH4_V(10))/(THKSL(10)*0.01)       !MODIFIED ON 2017 inserted  switch depend on hypothesis: the bottom boundary is a no-flux boundary or the 11th layer concentration is 0
        ! below I try to keep the CH4 flux no larger than the amount of CH4 that exist at the moment   V1.1 V1.2
        do i=1,nlayers+1
            if (Fdifu(i) .gt. 0.0 .and. (Fdifu(i)) .gt. CH4(i)) then
                Fdifu(i)=CH4(i)
            else if (Fdifu(i) .lt. 0.0 .and. (abs(Fdifu(i))) .gt. CH4(i-1)) then
                Fdifu(i)=-CH4(i-1)    
            endif
        enddo
          
        ! CH4(1) = CH4(1) + (0.0+Fdifu(1))/(THKSL(1)*0.01)
        do i = 1,nlayers-1                                  !loop of time
            CH4(i) = CH4(i) + (Fdifu(i+1)-Fdifu(i))*1 ! *1   * 1 hour /hour   /15min  *0.25h                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            ! CH4(i) = CH4(i) - 0.1*CH4(i)
            if (CH4(i) .lt. 0.0) then                     ! this part need to be improved until deleted   V1.2
                CH4(i) = 0.0
            endif
        enddo    
        CH4(10) = CH4(10) - Fdifu(10)                                   !MODIFIED ON 07/25/2016
        if (CH4(10) .lt. 0.0) then                                    !defined the Fdifu(11) to be 0.0
            CH4(10)= 0.0                                              ! switch on/off
        endif
        ! CH4(10) = CH4(10) +(Fdifu(11)-Fdifu(10))*1                      !MODIFIED ON 05/04/2018
        ! if (CH4(10) .lt. 0.0) then                                    
        !     CH4(10)= 0.0                                              
        ! endif        
        simuCH4 = simuCH4 + (Fdifu(1)-0.0) 
    
        ! ********************************************************************************************************************      
        ! D. methane ebullition     !assume bubbles can reach the water table within 1 h&
                                    !& the bubbles is added to the methane concentration in the soil layer just above the wt
                                    !& and then diffused through layers   ??not correct
        ! this subroutine is modified on 02132017 by deleting the unsat from bubble and add unsat to concentration so as to increase diffusion
        ! just by searching "switch" you can switch from old to new mode by adding or deleting "!"
        ! modified threshold value to 100 for testing
        ! ********************************************************************************************************************
        Kebu=1.0                    !unit  h-1   rate constant               
        Ebu_sum_unsat=0.0
        Ebu_sum_sat=0.0                                      !initial value
          
        do i=1,nlayers
            CH4_thre=1000.0  !!find in parafile  !unit  umol L-1 according to Walter's 500-1000
            CH4_thre_ly(i)=(CH4_thre*1.0e-6)*12*1000*(wsc(i)*0.001)    !convert the unit of CH4_thre from µmol L-1 to gC m-2
        enddo
        
        if (zwt .ge. 0.0) then                                  !when water table is above the soil surface
            do i=1,nlayers
                if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                    EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))     !only if the concentration is larger than threshold
                else !if (CH4(i) .le. CH4_thre_ly(i)) then
                    EbuCH4(i)=0.0
                endif
                Ebu_sum_sat=Ebu_sum_sat+EbuCH4(i)               !& the bubbles are directly added into CH4 efflux into atmosphere
                CH4(i)=CH4(i)- EbuCH4(i)                        !& update the concentration at the end of this hour in each layers
            enddo
        endif
        ! write (*,*) CH4(1),CH4_thre_ly(1),EbuCH4(1),Ebu_sum_sat 
        if (zwt .lt. 0.0) then                                  !when water table is below the soil surface
            do i=1,nlayers
                if ((depth(i)*10.0) .le. -zwt) then               !acrotelm layers
                    EbuCH4(i)=0.0
                    Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)         
                    CH4(i)=CH4(i)- EbuCH4(i) 
                else
                    if (((depth(i)*10.0)-(THKSL(i)*10.0)) .lt. -zwt) then       !partly acrotelm layer
                        wtlevelindex = i
                        if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                            EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))!*(((depth(i)*10.0)-(-zwt))/(THKSL(i)*10.0))        ! * percent
                        else !if (CH4(i) .le. CH4_thre_ly(i)) then                     ??????????,??????????????????
                            EbuCH4(i)=0.0
                        endif 
                        CH4(i)=CH4(i)- EbuCH4(i)
                        Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)                ! !  modified by Mary on 02132017
                        CH4(wtlevelindex-1)=CH4(wtlevelindex-1)+EbuCH4(i)    !!!!!-1-!!!! !switch on in new mode should be added add burst bubbles below surface to diffusion modified by Mary on 02132017
                        ! ：the problem is the resolution of soil layer is 10cm and EbuCH4(i) is directly added to the upper layer of boundary layer   02152017                    
                    else if (((depth(i)*10.0)-(THKSL(i)*10.0)) .ge. -zwt) then   !catotelm layers
                        if (CH4(i) .gt. CH4_thre_ly(i)) then                  
                            EbuCH4(i)=Kebu*(CH4(i)-CH4_thre_ly(i))
                        else !if (CH4(i) .le. CH4_thre_ly(i)) then
                            EbuCH4(i)=0.0
                        endif 
                        CH4(i)=CH4(i)- EbuCH4(i)        
                        CH4(wtlevelindex-1)=CH4(wtlevelindex-1)+EbuCH4(i)     !!!!!-2-!!!! !switch on in new mode should be added     modified by Mary on 02152017
                        Ebu_sum_unsat=Ebu_sum_unsat+EbuCH4(i)                  ! modified by Mary on 02132017
                    endif
                endif
            enddo
        endif
          
        Ebu_sum= Ebu_sum_sat
        ! simuCH4=simuCH4+Ebu_sum_sat                         !& the bubbles are directly added into CH4 efflux into atmosphere
        ! write (*,*) Ebu_sum        
        ! ******************************************************************************************************
        ! E. plant mediated methane transportation      totoally used Walter's model also used by Zhuang et. al
        ! ******************************************************************************************************
        Kpla=0.01         !unit h-1
        ! Kpla=0.01         !unit h-1
        ! Tveg=0.3 ! a factor describing the quality of plant-mediated transport depend on the density of plant stands and plant types 
        ! 0 for boreal forest and 0.5 for tundra
        ! find in parafile !
        ! the Tsoil used here would be better if refer to the 20cm soil temperature after &
        ! & the accomplishment of soil heat dynamics module. according to Zhuang. however Walter used 50cm soil temp.
        Tgr=2.0               !unit degree Celsius if annual mean temp is below 5 (otherwise 7)
        Tmat=Tgr+10.0         !unit degree Celsius
        Pox=0.5               !50% of mediated methane are oxidised 
        ! define fgrow
        if (Tsoil .lt. Tgr) then
            fgrow=LAIMIN
        else if (Tsoil .ge. Tgr .and. Tsoil .le. Tmat) then
            fgrow=LAIMIN+LAIMAX*(1-((Tmat-Tsoil)/(Tmat-Tgr))**2)
        else if (Tsoil .gt. Tmat) then
            fgrow=LAIMAX
        endif
          
        Pla_sum=0.0
        do i=1,nlayers
            PlaCH4(i)=Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)*(1-Pox)         !not sensitive at all to this change, but better
            ! PlaCH4(i)=Kpla*Tveg*FRLEN(i)*fgrow*CH4(i)
            Pla_sum=Pla_sum+PlaCH4(i)
            CH4(i)=CH4(i)-PlaCH4(i)
            CH4_V(i) = CH4(i)/(wsc(i)*0.001)          !convert concentration from gC/m2 to gC/m3
        enddo  
        simuCH4=simuCH4+Pla_sum
        consum=simuCH4+OxiCH4(1)+OxiCH4(2)+OxiCH4(3)+OxiCH4(4)+OxiCH4(5)+OxiCH4(6)+OxiCH4(7)+OxiCH4(8)+OxiCH4(9)+OxiCH4(10)
          
        !      if (MEMCMC .eq. 0) then
        !      !write(*,*) 'zwt',zwt,'simuCH4',simuCH4         !show on screen  
        ! ***********     write out hourly value for methane module 
    
        !      write(82,182)zwt,Pla_sum,simuCH4, &
        !              & Rh(1),Rh(2),Rh(3),Rh(4),Rh(5),Rh(6),Rh(7),Rh(8),Rh(9),Rh(10),   &
        !              & consum,Pro_sum, &
        !              & ProCH4(1),ProCH4(2),ProCH4(3),ProCH4(4),ProCH4(5),ProCH4(6),ProCH4(7),ProCH4(8),ProCH4(9),ProCH4(10),   &
        !              & CH4(1),CH4(2),CH4(3),CH4(4),CH4(5),CH4(6),CH4(7),CH4(8),CH4(9),CH4(10), &
        !              & CH4_V(1),CH4_V(2),CH4_V(3),CH4_V(4),CH4_V(5),CH4_V(6),CH4_V(7),CH4_V(8),CH4_V(9),CH4_V(10), &              
        !              & Fdifu(1),Fdifu(2),Fdifu(3),Fdifu(4),Fdifu(5),Fdifu(6),Fdifu(7),Fdifu(8),Fdifu(9),Fdifu(10),Fdifu(11), &
        !              & OxiCH4(1),OxiCH4(2),OxiCH4(3),OxiCH4(4),OxiCH4(5),OxiCH4(6),OxiCH4(7),OxiCH4(8),OxiCH4(9),OxiCH4(10),   &
        !              & wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10), &
        !              & Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Tsoil
        !
        !182   format(81(f15.9,","))  
        !          Ebu_sum_unsat=0.0
        !          Ebu_sum_sat=0.0 
        !           write (*,*) Ebu_sum_sat,Ebu_sum_unsat
        !        write(82,182) zwt, Rh(1), Rh_pools(1)
        !182     format(5(f15.9,","))  
        !        write (121,1201) Rh_pools(1) !Ebu_sum_sat! Ebu_sum_unsat
        !1201    format(2(f15.4,","))        
        !            consum,Pro_sum,zwt,Ebu_sum_sat,Rh(1),Rh(2),Rh(3),Rh(4),Rh(5),Rh(6),Rh(7),Rh(8),Rh(9),Rh(10),  &
        !        &   ProCH4(1),ProCH4(2),ProCH4(3),ProCH4(4),ProCH4(5),ProCH4(6),ProCH4(7),ProCH4(8),ProCH4(9),ProCH4(10), &
        !        &   CH4(1),CH4(2),CH4(3),CH4(4),CH4(5),CH4(6),CH4(7),CH4(8),CH4(9),CH4(10),  &
        !        &   ProCH4(1),ProCH4(2),ProCH4(3),ProCH4(4),ProCH4(5),ProCH4(6),ProCH4(7),ProCH4(8),ProCH4(9),ProCH4(10), &
        !        &   Fdifu(1),Fdifu(2),Fdifu(3),Fdifu(4),Fdifu(5),Fdifu(6),Fdifu(7),Fdifu(8),Fdifu(9),Fdifu(10),  &
        !        &   OxiCH4(1),OxiCH4(2),OxiCH4(3),OxiCH4(4),OxiCH4(5),OxiCH4(6),OxiCH4(7),OxiCH4(8),OxiCH4(9),OxiCH4(10),  &
        !        &   wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10),  &
        !        &   Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Tsoil
        !      write(83,183)zwt,Tsoil,Rh_pools(1),Rh_pools(2),Rh_pools(3),Rh_pools(4),Rh_pools(5), &
        !              & wsc(1),wsc(2),wsc(3),wsc(4),wsc(5),wsc(6),wsc(7),wsc(8),wsc(9),wsc(10)
        !
        !183   format(17(f15.9,","))      
        !!      endif    
        ! 
        ! ***********     write out hourly value for methane module 
        return
    end subroutine methane

    real function esat(T)   ! returns saturation vapour pressure in Pa
        real T
        esat = 610.78*exp(17.27*T/(T+237.3))
        return
    end

end module mod_soil