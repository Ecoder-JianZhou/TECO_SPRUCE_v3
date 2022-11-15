module mod_teco
    implicit none

    contains
    subroutine canopy()
        ! Normalised Gaussian points and weights (Goudriaan & van Laar, 1993, P98)
        ! 5-point
        data Gaussx/0.0469101,0.2307534,0.5,0.7692465,0.9530899/
        data Gaussw/0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/
        data Gaussw_cum/0.11846,0.35777,0.64222,0.88153,1.0/
        ! calculate beam fraction in incoming solar radiation
        call  yrday(doy,hour,lat,radsol,fbeam)
        idoy=int(doy)
        hours=idoy*1.0+hour/24.0
        coszen=sinbet(doy,lat,pi,hour)             !cos zenith angle of sun
        ! set windspeed to the minimum speed to avoid zero Gb
        if(wind.lt.0.01) wind=0.01
        ! calculate soil albedo for NIR as a function of soil water (Garratt pp292)
        if(topfws.gt.0.5) then
            rhoS(2)=0.18
        else
            rhoS(2)=0.52-0.68*topfws
        endif
        ! assign plant biomass and leaf area index at time t
        ! assume leaf biomass = root biomass
        FLAIT =LAI 
        eairP=esat(Tair)-Dair                !air water vapour pressure
        radabv(1)=0.5*radsol                 !(1) - solar radn
        radabv(2)=0.5*radsol                 !(2) - NIR
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
end module mod_teco