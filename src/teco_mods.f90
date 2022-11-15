module mod_teco
    ! each loop for simulation
    use mod_data
    implicit none
    integer doy, hour
    real radsol, wind, tair, VPD, TairK, co2ca
    real coszen, fbeam, flait, Radabv(2), eairP
    real Acan1, Acan2, Ecan1, Ecan2
    real raero, extKb, extkd
    real flai, Qabs(3,2), emair, Rnstar(2), grdn
    real reff(3,2),kpr(3,2),scatt(2)
    real windUx, Vcmxx, eJmxx
    real Gaussx(5),Gaussw(5),Gaussw_cum(5)                          ! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
    data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/        ! 5-point
    data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
    data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/

    contains
    subroutine canopy(in_doy,in_hour,in_radsol, in_wind, in_tair, in_VPD, in_co2ca)
        implicit none
        integer in_doy, in_hour
        real in_radsol, in_wind, in_tair, in_VPD, in_co2ca
        
        doy    = in_doy
        hour   = in_hour
        radsol = in_radsol
        wind   = in_wind
        tair   = in_tair
        VPD    = in_VPD             ! Dair
        co2ca  = in_co2ca           ! co2ca = 380.0*1.0E-6
        
        call  yrday()              ! calculate beam fraction in incoming solar radiation
        ! hours  = int(doy)*1.0+hour/24.0               ! Jian: seem no used
        coszen = sinbet()             !cos zenith angle of sun
        ! set windspeed to the minimum speed to avoid zero Gb
        if(wind.lt.0.01) wind=0.01
        ! calculate soil albedo for NIR as a function of soil water (Garratt pp292)
        if(topfws.gt.0.5) then
            rhoS(2) = 0.18
        else
            rhoS(2) = 0.52-0.68*topfws
        endif
        ! assign plant biomass and leaf area index at time t
        ! assume leaf biomass = root biomass
        FLAIT     = LAI 
        eairP     = esat(Tair)-VPD                !air water vapour pressure
        radabv(1) = 0.5*radsol                 !(1) - solar radn
        radabv(2) = 0.5*radsol                 !(2) - NIR
        ! call multilayer model of Leuning - uses Gaussian integration but radiation scheme
        ! is that of Goudriaan
        call xlayers(Sps,Tair,Dair,radabv,fbeam,eairP,&                                   ! 
            &        wind,co2ca,fwsoil,wcl,FLAIT,coszen,idoy,hours,&
            &        tauL,rhoL,rhoS,xfang,extkd,extkU,wleaf,&
            &        Rconst,sigma,emleaf,emsoil,theta,a1,Ds0,&
            &        cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,&
            &        gsw0,alpha,stom_n,wsmax,wsmin,&
            &        Vcmx0,eJmx0,conKc0,conKo0,Ekc,Eko,o2ci,&
            &        Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,&
            &        extKb,Rsoilabs,Acan1,Acan2,Ecan1,Ecan2,&
            &        RnStL,QcanL,RcanL,AcanL,EcanL,HcanL,GbwcL,GswcL,gddonset,&
            &        testout,Rsoilab1,Rsoilab2,QLleaf,QLair,raero,do_soilphy,&  ! added from soil thermal ..int 
            &        G,Esoil,Hsoil) ! added from soil thermal ..int 
        ! *** ..int   added 'testout,Rsoilab1,Rsoilab2,QLleaf,QLair,raero,do_soilphy,G,Esoil,Hsoil'
        if (do_soilphy) then 
            call Tsoil_simu(Rsoilab1,Rsoilab2,QLleaf,QLair,Tair,Dair,&
                &         fbeam,FLAIT,sigma,emsoil,rhoS,Rconst,&
                &         extkd,extkb,cpair,Patm,AirMa,H2OMw,&
                &         H2OLv0,wcl,raero,wsmax,wsmin,wind,sftmp,Tsoill,testout,ht,ice,&
                &         snow_depth,Tsnow,Twater,Tice,water_tw,ice_tw,diff_s,G,tsoil,&
                &         diff_snow,albedo_snow,resht,thd_snow_depth,thksl,zwt,Esoil,Hsoil,liq_water,shcap_snow,&
                &         condu_snow,condu_b,depth_ex,dcount_soil)
        endif
        !     write (84,184) Esoil
        !184   format(f15.9,",")
        !   ***  
        Acanop=Acan1+Acan2
        Ecanop=Ecan1+Ecan2
        gpp=Acanop*3600.0*12.0                           ! every hour, g C m-2 h-1
        transp=AMAX1(Ecanop*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.) ! mm H2O /hour
        evap=AMAX1(Esoil*3600.0/(1.0e6*(2.501-0.00236*Tair)),0.)
        ! evap=evap*0.8
        ! H2OLv0=2.501e6               !latent heat H2O (J/kg)
        return
    end subroutine canopy

    subroutine yrday()
        real pidiv, slatx, sindec, cosdec
        real a, b, sinbet0, solext, tmprat, tmpR, tmpK, fdiff
        pidiv  = pi/180.0
        slatx  = lat*pidiv
        sindec = -sin(23.4*pidiv)*cos(2.0*pi*(doy+10.0)/365.0)
        cosdec = sqrt(1.-sindec*sindec)
        a      = sin(slatx)*sindec
        b      = cos(slatx)*cosdec
        sinbet0 = a+b*cos(2*pi*(hour-12.)/24.)
        solext = 1370.0*(1.0+0.033*cos(2.0*pi*(doy-10.)/365.0))*sinbet0
        tmprat = radsol/solext
        tmpR   = 0.847-1.61*sinbet0+1.04*sinbet0*sinbet0
        tmpK   = (1.47-tmpR)/1.66
        if(tmprat.le.0.22) fdiff=1.0
        if(tmprat.gt.0.22.and.tmprat.le.0.35) then
            fdiff = 1.0-6.4*(tmprat-0.22)*(tmprat-0.22)
        endif
        if(tmprat.gt.0.35.and.tmprat.le.tmpK) then
            fdiff = 1.47-1.66*tmprat
        endif
        if(tmprat.ge.tmpK) then
            fdiff = tmpR
        endif
        fbeam = 1.0-fdiff
        if(fbeam.lt.0.0) fbeam = 0.0
        return
    end subroutine yrday

! used by canopy module ------------------------------------------------------
    subroutine xlayers() ! added from soil thermal ..int 
        ! the multi-layered canopy model developed by 
        ! Ray Leuning with the new radiative transfer scheme   
        ! implemented by Y.P. Wang (from Sellers 1986)
        ! 12/Sept/96 (YPW) correction for mean surface temperature of sunlit
        ! and shaded leaves
        ! Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)} 
        ! ---------------------------------------------------------------------
        ! real layer1(5),layer2(5)
        ! real tauL(3),rhoL(3),rhoS(3),Qabs(3,2),Radabv(2),Rnstar(2)
        ! real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
        ! real gbleaf(2),gsleaf(2),QSabs(3,2),Qasoil(2)
        ! integer ng,nw
        ! real rhoc(3,2),reff(3,2),kpr(3,2),scatt(2)       !Goudriaan
    
        ! real rsoil,rlai,raero,LAI
        ! real wsmax,wsmin,WILTPT,FILDCP,wcl(10)
        ! real gddonset
        ! ! additional arrays to allow output of info for each Layer
        ! real RnStL(5),QcanL(5),RcanL(5),AcanL(5),EcanL(5),HcanL(5)
        ! real GbwcL(5),GswcL(5)
    
        ! !   *** ..int
        ! !*************************
        ! real testout(11)      
        ! logical do_soilphy
        !   *** .int
        real :: Rnst1 = 0.0,  Rnst2 = 0.0, Qcan1 = 0.0, Qcan2 = 0.0 !net rad, sunlit, vis rad
        real :: Rcan1 = 0.0,  Rcan2 = 0.0, Hcan1 = 0.0, Hcan2 = 0.0 !NIR rad !Sens heat
        real :: Gbwc1 = 0.0,  Gbwc2 = 0.0, Gswc1 = 0.0, Gswc2 = 0.0 !Boundary layer conductance; Canopy conductance
        real :: Tleaf1 = 0.0, Tleaf2 = 0.0                          !Leaf Temp
        real xphi1, xphi2, funG ! more places?
        real pi180, cozen15, cozen45, cozen75, xK15, xK45, xK75
        real transd, extkn
        integer nw, ng
        real rhoc(3,2)       !Goudriaan
        real rhoch, rhoc15, rhoc45, rhoc75
        real scalex
        ! soil water conditions
        WILTPT=wsmin/100.
        FILDCP=wsmax/100.
        ! reset the vairables
        Acan1=0.0        !CO2
        Acan2=0.0
        Ecan1=0.0        !Evap
        Ecan2=0.0
        ! aerodynamic resistance                                                
        raero=50./wind                           
        ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
        xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
        xphi2 = 0.877 * (1.0 - 2.0*xphi1)
        funG  = xphi1 + xphi2*coszen                             !G-function: Projection of unit leaf area in direction of beam
        if(coszen.gt.0) then                                     !check if day or night
            extKb=funG/coszen                                    !beam extinction coeff - black leaves
        else
            extKb=100.
        end if
        ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
        ! Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
        pi180   = 3.1416/180.
        cozen15 = cos(pi180*15)
        cozen45 = cos(pi180*45)
        cozen75 = cos(pi180*75)
        xK15    = xphi1/cozen15+xphi2
        xK45    = xphi1/cozen45+xphi2
        xK75    = xphi1/cozen75+xphi2
        transd  = 0.308*exp(-xK15*FLAIT)+0.514*exp(-xK45*FLAIT)+     &
                  &       0.178*exp(-xK75*FLAIT)
        extkd   = (-1./FLAIT)*alog(transd)
        extkn   = extkd                        !N distribution coeff 
    
        ! canopy reflection coefficients (Array indices: first;  1=VIS,  2=NIR
                                                        !   second; 1=beam, 2=diffuse
        do nw=1,2  ! nw:1=VIS, 2=NIR
            scatt(nw)  = tauL(nw)+rhoL(nw)                                          ! scattering coeff
            if((1.-scatt(nw))<0.0)scatt(nw)=0.9999                                  ! Weng 10/31/2008
            kpr(nw,1)  = extKb*sqrt(1.-scatt(nw))                                   ! modified k beam scattered (6.20)
            kpr(nw,2)  = extkd*sqrt(1.-scatt(nw))                                   ! modified k diffuse (6.20)
            rhoch      = (1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))            ! canopy reflection black horizontal leaves (6.19)
            rhoc15     = 2.*xK15*rhoch/(xK15+extkd)                                 ! canopy reflection (6.21) diffuse
            rhoc45     = 2.*xK45*rhoch/(xK45+extkd)
            rhoc75     = 2.*xK75*rhoch/(xK75+extkd)
            rhoc(nw,2) = 0.308*rhoc15+0.514*rhoc45+0.178*rhoc75
            rhoc(nw,1) = 2.*extKb/(extKb+extkd)*rhoch                               ! canopy reflection (6.21) beam 
            reff(nw,1) = rhoc(nw,1)+(rhoS(nw)-rhoc(nw,1))*exp(-2.*kpr(nw,1)*FLAIT)  ! effective canopy-soil reflection coeff - beam (6.27)          
            reff(nw,2) = rhoc(nw,2)+(rhoS(nw)-rhoc(nw,2))*exp(-2.*kpr(nw,2)*FLAIT)  ! effective canopy-soil reflection coeff - diffuse (6.27)         
        enddo
        ! isothermal net radiation & radiation conductance at canopy top - needed to calc emair
        call Radiso()           ! Jian: some parameters not initialization.
        TairK = Tair+273.2
        ! below      
        do ng=1,5
            flai = gaussx(ng)*FLAIT
            ! radiation absorption for visible and near infra-red
            call goudriaan() 
            ! isothermal net radiation & radiation conductance at canopy top
            call Radiso(flai,flait,Qabs,extkd,Tair,eairP,cpair,Patm,   &
                &               fbeam,airMa,Rconst,sigma,emleaf,emsoil,        &
                &               emair,Rnstar,grdn)
            windUx = wind*exp(-extkU*flai)                ! windspeed at depth xi
            scalex = exp(-extkn*flai)                     ! scale Vcmx0 & Jmax0
            Vcmxx  = Vcmax0*scalex                         ! Vcmx0 ---> Vcmax0
            eJmxx  = eJmx0*scalex
            if(radabv(1).ge.10.0) then                          !check solar Radiation > 10 W/m2
            ! leaf stomata-photosynthesis-transpiration model - daytime
                call agsean_day()       
            else
                call agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,      &
                    &    co2ca,wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,  &
                    &    Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,  &
                    &    gsw0,alpha,stom_n,                                 &
                    &    Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,            &
                    &    Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,         &
                    &    Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf)
            endif  
            fslt=exp(-extKb*flai)                        !fraction of sunlit leaves
            fshd=1.0-fslt                                !fraction of shaded leaves
            Rnst1=Rnst1+fslt*Rnstar(1)*Gaussw(ng)*FLAIT  !Isothermal net rad`
            Rnst2=Rnst2+fshd*Rnstar(2)*Gaussw(ng)*FLAIT
            RnstL(ng)=Rnst1+Rnst2
    
            Qcan1=Qcan1+fslt*Qabs(1,1)*Gaussw(ng)*FLAIT  !visible
            Qcan2=Qcan2+fshd*Qabs(1,2)*Gaussw(ng)*FLAIT
            QcanL(ng)=Qcan1+Qcan2
    
            Rcan1=Rcan1+fslt*Qabs(2,1)*Gaussw(ng)*FLAIT  !NIR
            Rcan2=Rcan2+fshd*Qabs(2,2)*Gaussw(ng)*FLAIT
            RcanL(ng)=Rcan1+Rcan2
    
            if(Aleaf(1).lt.0.0)Aleaf(1)=0.0      !Weng 2/16/2006
            if(Aleaf(2).lt.0.0)Aleaf(2)=0.0      !Weng 2/16/2006
    
            Acan1=Acan1+fslt*Aleaf(1)*Gaussw(ng)*FLAIT*stom_n    !amphi/hypostomatous
            Acan2=Acan2+fshd*Aleaf(2)*Gaussw(ng)*FLAIT*stom_n
            AcanL(ng)=Acan1+Acan2
    
            layer1(ng)=Aleaf(1)
            layer2(ng)=Aleaf(2)
    
            Ecan1=Ecan1+fslt*Eleaf(1)*Gaussw(ng)*FLAIT
            Ecan2=Ecan2+fshd*Eleaf(2)*Gaussw(ng)*FLAIT
            EcanL(ng)=Ecan1+Ecan2
    
            Hcan1=Hcan1+fslt*Hleaf(1)*Gaussw(ng)*FLAIT
            Hcan2=Hcan2+fshd*Hleaf(2)*Gaussw(ng)*FLAIT
            HcanL(ng)=Hcan1+Hcan2
    
            Gbwc1=Gbwc1+fslt*gbleaf(1)*Gaussw(ng)*FLAIT*stom_n
            Gbwc2=Gbwc2+fshd*gbleaf(2)*Gaussw(ng)*FLAIT*stom_n
    
            Gswc1=Gswc1+fslt*gsleaf(1)*Gaussw(ng)*FLAIT*stom_n
            Gswc2=Gswc2+fshd*gsleaf(2)*Gaussw(ng)*FLAIT*stom_n
    
            Tleaf1=Tleaf1+fslt*Tleaf(1)*Gaussw(ng)*FLAIT
            Tleaf2=Tleaf2+fshd*Tleaf(2)*Gaussw(ng)*FLAIT
        enddo  ! 5 layers
    
        FLAIT1=(1.0-exp(-extKb*FLAIT))/extkb
        Tleaf1=Tleaf1/FLAIT1
        Tleaf2=Tleaf2/(FLAIT-FLAIT1)
        ! Soil surface energy and water fluxes
        ! Radiation absorbed by soil
        Rsoilab1=fbeam*(1.-reff(1,1))*exp(-kpr(1,1)*FLAIT)        &
            &         +(1.-fbeam)*(1.-reff(1,2))*exp(-kpr(1,2)*FLAIT)          !visible
        Rsoilab2=fbeam*(1.-reff(2,1))*exp(-kpr(2,1)*FLAIT)        &
            &         +(1.-fbeam)*(1.-reff(2,2))*exp(-kpr(2,2)*FLAIT)          !NIR
        Rsoilab1=Rsoilab1*Radabv(1)
        Rsoilab2=Rsoilab2*Radabv(2)  
        Tlk1=Tleaf1+273.2
        Tlk2=Tleaf2+273.2
        ! temp1=-extkd*FLAIT
        QLair=emair*sigma*(TairK**4)*exp(-extkd*FLAIT)
        QLleaf=emleaf*sigma*(Tlk1**4)*exp(-extkb*FLAIT)           &
            &      +emleaf*sigma*(Tlk2**4)*(1.0-exp(-extkb*FLAIT))
        QLleaf=QLleaf*(1.0-exp(-extkd*FLAIT)) 
        QLsoil=emsoil*sigma*(TairK**4)
        Rsoilab3=(QLair+QLleaf)*(1.0-rhoS(3))-QLsoil
    
        ! Net radiation absorbed by soil
        ! the old version of net long-wave radiation absorbed by soils 
        ! (with isothermal assumption)
        ! Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*FLAIT)         !Longwave
        ! Rsoilab3=(1-rhoS(3))*Rsoil3
    
        ! Total radiation absorbed by soil    
        Rsoilabs=Rsoilab1+Rsoilab2+Rsoilab3 
    
        ! thermodynamic parameters for air
        TairK=Tair+273.2
        rhocp=cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv=H2oLv0-2.365e3*Tair
        slope=(esat(Tair+0.1)-esat(Tair))/0.1
        psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar=Patm/(Rconst*TairK)
        fw1=AMIN1(AMAX1((FILDCP-wcl(1))/(FILDCP-WILTPT),0.05),1.0)
        Rsoil=30.*exp(0.2/fw1)
        rLAI=exp(FLAIT)
        ! latent heat flux into air from soil
        ! Eleaf(ileaf)=1.0*
        ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
        ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
        Esoil=(slope*(Rsoilabs-G)+rhocp*Dair/(raero+rLAI))/       &
            &      (slope+psyc*(rsoil/(raero+rLAI)+1.))
        ! sensible heat flux into air from soil
        Hsoil=Rsoilabs-Esoil-G
        return
    end 

    subroutine Radiso()
        ! output
        ! Rnstar(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
        ! 23 Dec 1994
        ! calculates isothermal net radiation for sunlit and shaded leaves under clear skies
        ! implicit real (a-z)
        real Rnstar(2)
        real Qabs(3,2)
        TairK=Tair+273.2
    
        ! thermodynamic properties of air
        rhocp=cpair*Patm*airMa/(Rconst*TairK)   !volumetric heat capacity (J/m3/K)
    
        ! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
        emsky=0.642*(eairP/Tairk)**(1./7)       !note eair in Pa
         
        ! apparent emissivity from clouds (Kimball et al 1982)
        ep8z=0.24+2.98e-12*eairP*eairP*exp(3000/TairK)
        tau8=amin1(1.0,1.0-ep8z*(1.4-0.4*ep8z))            !ensure tau8<1
        emcloud=0.36*tau8*(1.-fbeam)*(1-10./TairK)**4      !10 from Tcloud = Tair-10
    
        ! apparent emissivity from sky plus clouds      
        !      emair=emsky+emcloud
        ! 20/06/96
        emair=emsky
        if(emair.gt.1.0) emair=1.0
        ! net isothermal outgoing longwave radiation per unit leaf area at canopy
        ! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
        Bn0=sigma*(TairK**4.)
        Bnxi=Bn0*extkd*(exp(-extkd*flai)*(emair-emleaf)       &
            &    + exp(-extkd*(flait-flai))*(emsoil-emleaf))
        ! isothermal net radiation per unit leaf area for thin layer of sunlit and
        ! shaded leaves
        Rnstar(1)=Qabs(1,1)+Qabs(2,1)+Bnxi
        Rnstar(2)=Qabs(1,2)+Qabs(2,2)+Bnxi
        ! radiation conductance (m/s) @ flai
        grdn=4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
            &    (exp(-extkd*flai)+exp(-extkd*(flait-flai)))       &
            &    /rhocp
        return
    end

    subroutine goudriaan()
        ! for spheric leaf angle distribution only
        ! compute within canopy radiation (PAR and near infra-red bands)
        ! using two-stream approximation (Goudriaan & vanLaar 1994)
        ! tauL: leaf transmittance
        ! rhoL: leaf reflectance
        ! rhoS: soil reflectance
        ! sfang XiL function of Ross (1975) - allows for departure from spherical LAD
        !     (-1 vertical, +1 horizontal leaves, 0 spherical)
        ! FLAI: canopy leaf area index
        ! funG: Ross' G function
        ! scatB: upscatter parameter for direct beam
        ! scatD: upscatter parameter for diffuse
        ! albedo: single scattering albedo
        ! output:
        ! Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
        !                     type=1 for sunlit;   =2 for shaded (W/m2)
        real xu, xphi1, xphi2, funG, Qd0, Qb0
        integer nw
        xu = coszen                                         !cos zenith angle
          
        ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
        xphi1 = 0.5 - 0.633*xfang -0.33*xfang*xfang
        xphi2 = 0.877 * (1.0 - 2.0*xphi1)
        funG  = xphi1 + xphi2*xu                             !G-function: Projection of unit leaf area in direction of beam
        if(coszen.gt.0) then                                  !check if day or night
            extKb=funG/coszen                                   !beam extinction coeff - black leaves
        else
            extKb=100.
        end if                 
        ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
        do nw=1,2
            Qd0=(1.-fbeam)*radabv(nw)                                          !diffuse incident radiation
            Qb0=fbeam*radabv(nw)                                               !beam incident radiation
            Qabs(nw,2)=Qd0*(kpr(nw,2)*(1.-reff(nw,2))*exp(-kpr(nw,2)*FLAI))+  & !absorbed radiation - shaded leaves, diffuse
                &      Qb0*(kpr(nw,1)*(1.-reff(nw,1))*exp(-kpr(nw,1)*FLAI)-   & !beam scattered
                &      extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
            Qabs(nw,1)=Qabs(nw,2)+extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves 
        end do
        return
    end

    subroutine agsean_day()
        integer kr1,ileaf
        real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
        real gbleaf(2), gsleaf(2)
        real Qabs(3,2),Rnstar(2)
        ! thermodynamic parameters for air
        TairK=Tair+273.2
        rhocp=cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv=H2oLv0-2.365e3*Tair
        slope=(esat(Tair+0.1)-esat(Tair))/0.1
        psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar=Patm/(Rconst*TairK)
        weighJ=1.0
        ! boundary layer conductance for heat - single sided, forced convection
        ! (Monteith 1973, P106 & notes dated 23/12/94)
        if(windUx/wleaf>=0.0)then
            gbHu=0.003*sqrt(windUx/wleaf)    !m/s
        else
            gbHu=0.003 !*sqrt(-windUx/wleaf)
        endif         ! Weng 10/31/2008
        ! raero=0.0                        !aerodynamic resistance s/m
        do ileaf=1,2              ! loop over sunlit and shaded leaves
            ! first estimate of leaf temperature - assume air temp
            Tleaf(ileaf)=Tair
            Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
            ! first estimate of deficit at leaf surface - assume Da
            Dleaf=Dair                !Pa
            ! first estimate for co2cs
            co2cs=co2ca               !mol/mol
            Qapar = (4.6e-6)*Qabs(1,ileaf)
            ! ********************************************************************
            kr1=0                     !iteration counter for LE
            ! return point for evaporation iteration
            do               !iteration for leaf temperature
                ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
                Gras=1.595e8*ABS(Tleaf(ileaf)-Tair)*(wleaf**3.)     !Grashof
                gbHf=0.5*Dheat*(Gras**0.25)/wleaf
                gbH=gbHu+gbHf                         !m/s
                rbH=1./gbH                            !b/l resistance to heat transfer
                rbw=0.93*rbH                          !b/l resistance to water vapour
                ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
                rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
                rrdn=1./grdn
                Y=1./(1.+ (rbH_L+raero)/rrdn)
                ! boundary layer conductance for CO2 - single side only (mol/m2/s)
                gbc=Cmolar*gbH/1.32            !mol/m2/s
                gsc0=gsw0/1.57                 !convert conductance for H2O to that for CO2
                varQc=0.0
                weighR=1.0
                call photosyn(Sps,CO2Ca,CO2Cs,Dleaf,Tlk,Qapar,Gbc,   &   !Qaparx<-Qapar,Gbcx<-Gsc0
                    &         theta,a1,Ds0,fwsoil,varQc,weighR,                &
                    &         gsc0,alpha,Vcmxx,eJmxx,weighJ,                   &
                    &         conKc0,conKo0,Ekc,Eko,o2ci,Rconst,Trefk,         &
                    &         Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,       &
                    &         Aleafx,Gscx,gddonset)  !outputs
                ! choose smaller of Ac, Aq
                Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
                ! calculate new values for gsc, cs (Lohammer model)
                co2cs = co2ca-Aleaf(ileaf)/gbc
                co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gscx
                ! scale variables
                ! gsw=gscx*1.56      !gsw in mol/m2/s, oreginal:gsw=gsc0*1.56,Weng20060215
                gsw=gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
                gswv=gsw/Cmolar                           !gsw in m/s
                rswv=1./gswv
                ! calculate evap'n using combination equation with current estimate of gsw
                Eleaf(ileaf)=1.0*(slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    &   !2* Weng 0215
                    &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))        
                ! calculate sensible heat flux
                Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
                ! calculate new leaf temperature (K)
                Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
                ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
                Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
                gbleaf(ileaf)=gbc*1.32*1.075
                gsleaf(ileaf)=gsw
                ! compare current and previous leaf temperatures
                if(abs(Tlk1-Tlk).le.0.1) exit ! original is 0.05 C Weng 10/31/2008
                ! update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
                Tlk=Tlk1
                Tleaf(ileaf)=Tlk1-273.2
                kr1=kr1+1
                if(kr1 > 500)then
                    Tlk=TairK
                    exit
                endif
                if(Tlk < 200.)then
                    Tlk=TairK
                    exit 
                endif                     ! Weng 10/31/2008
                ! goto 100                          !solution not found yet
            enddo
    ! 10  continue
        enddo
        return
    end
    
    ! ****************************************************************************
    subroutine agsean_ngt(Sps,Qabs,Rnstar,grdn,windUx,Tair,Dair,co2ca,    &
        &               wleaf,raero,theta,a1,Ds0,fwsoil,idoy,hours,            &
        &               Rconst,cpair,Patm,Trefk,H2OLv0,AirMa,H2OMw,Dheat,      &
        &               gsw0,alpha,stom_n,                                     &
        &               Vcmxx,eJmxx,conKc0,conKo0,Ekc,Eko,o2ci,                &
        &               Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2,             &
        &               Aleaf,Eleaf,Hleaf,Tleaf,gbleaf,gsleaf,co2ci)
        ! implicit real (a-z)
        integer kr1,ileaf
        real Aleaf(2),Eleaf(2),Hleaf(2),Tleaf(2),co2ci(2)
        real gbleaf(2), gsleaf(2)
        real Qabs(3,2),Rnstar(2)
        ! thermodynamic parameters for air
        TairK=Tair+273.2
        rhocp=cpair*Patm*AirMa/(Rconst*TairK)
        H2OLv=H2oLv0-2.365e3*Tair
        slope=(esat(Tair+0.1)-esat(Tair))/0.1
        psyc=Patm*cpair*AirMa/(H2OLv*H2OMw)
        Cmolar=Patm/(Rconst*TairK)
        weighJ=1.0
        ! boundary layer conductance for heat - single sided, forced convection
        ! (Monteith 1973, P106 & notes dated 23/12/94)
        gbHu=0.003*sqrt(windUx/wleaf)    !m/s
        ! raero=0.0                        !aerodynamic resistance s/m
        do ileaf=1,2                  ! loop over sunlit and shaded leaves
            ! first estimate of leaf temperature - assume air temp
            Tleaf(ileaf)=Tair
            Tlk=Tleaf(ileaf)+273.2    !Tleaf to deg K
            ! first estimate of deficit at leaf surface - assume Da
            Dleaf=Dair                !Pa
            ! first estimate for co2cs
            co2cs=co2ca               !mol/mol
            Qapar = (4.6e-6)*Qabs(1,ileaf)
            ! ********************************************************************
            kr1=0                     !iteration counter for LE
            do
    !100        continue !    return point for evaporation iteration
                ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
                Gras=1.595e8*abs(Tleaf(ileaf)-Tair)*(wleaf**3)     !Grashof
                gbHf=0.5*Dheat*(Gras**0.25)/wleaf
                gbH=gbHu+gbHf                         !m/s
                rbH=1./gbH                            !b/l resistance to heat transfer
                rbw=0.93*rbH                          !b/l resistance to water vapour
                ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
                rbH_L=rbH*stom_n/2.                   !final b/l resistance for heat  
                rrdn=1./grdn
                Y=1./(1.+ (rbH_L+raero)/rrdn)
                ! boundary layer conductance for CO2 - single side only (mol/m2/s)
                gbc=Cmolar*gbH/1.32            !mol/m2/s
                gsc0=gsw0/1.57                        !convert conductance for H2O to that for CO2
                varQc=0.0                  
                weighR=1.0
                ! respiration      
                Aleafx=-0.0089*Vcmxx*exp(0.069*(Tlk-293.2))
                gsc=gsc0
                ! choose smaller of Ac, Aq
                Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
                ! calculate new values for gsc, cs (Lohammer model)
                co2cs = co2ca-Aleaf(ileaf)/gbc
                co2Ci(ileaf) = co2cs-Aleaf(ileaf)/gsc
                ! scale variables
                gsw=gsc*1.56                              !gsw in mol/m2/s
                gswv=gsw/Cmolar                           !gsw in m/s
                rswv=1./gswv
                ! calculate evap'n using combination equation with current estimate of gsw
                Eleaf(ileaf)= (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/   &
                    &      (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
                ! calculate sensible heat flux
                Hleaf(ileaf)=Y*(Rnstar(ileaf)-Eleaf(ileaf))
                ! calculate new leaf temperature (K)
                Tlk1=273.2+Tair+Hleaf(ileaf)*(rbH/2.+raero)/rhocp
                ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
                Dleaf=psyc*Eleaf(ileaf)/(rhocp*gswv)
                gbleaf(ileaf)=gbc*1.32*1.075
                gsleaf(ileaf)=gsw
    
                ! compare current and previous leaf temperatures
                if(abs(Tlk1-Tlk).le.0.1)exit
                if(kr1.gt.500)exit
                ! update leaf temperature
                Tlk=Tlk1 
                Tleaf(ileaf)=Tlk1-273.2
                kr1=kr1+1
            enddo                          !solution not found yet
    10    continue
        enddo
        return
    end


!  functions used in canopy -------------------------------
    real function sinbet()
        real rad, sinlat, coslat, sindec, cosdec, A, B
        ! sin(bet), bet = elevation angle of sun
        ! calculations according to Goudriaan & van Laar 1994 P30
        rad = pi/180.
        ! sine and cosine of latitude
        sinlat = sin(rad*lat)
        coslat = cos(rad*lat)
        ! sine of maximum declination
        sindec=-sin(23.45*rad)*cos(2.0*pi*(doy+10.0)/365.0)
        cosdec=sqrt(1.-sindec*sindec)
        ! terms A & B in Eq 3.3
        A = sinlat*sindec
        B = coslat*cosdec
        sinbet = A+B*cos(pi*(hour-12.)/12.)
        return
    end

    real function esat(T)
        real T
        ! returns saturation vapour pressure in Pa
        esat=610.78*exp(17.27*T/(T+237.3))
        return
    end
end module mod_teco