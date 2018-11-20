subroutine ncwrt_onlysed(iogrp)

!-----------------------------------------------------------------------
!
! Write diagnostic sediment fields during off-line sediment spin-up
!
! Copyright (C) 2018 Marco van Hulten <Marco.Hulten@uib.no> et al.
!                    Geophysical Institute @ University of Bergen
!
! This subroutine is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!
! Method:
!
! History:
!
! This is based on ncout_hamocc.F that contains ncwrt_bgc().
!
! Marco van Hulten <Marco.Hulten@uib.no>           2018-11-20
! - Initial code based on ncout_hamocc.F.
!
!-----------------------------------------------------------------------

use mod_xc
use mod_dia       , only: diafnm,sigmar1,iotype
use mo_control_bgc, only: dtbgc, nburst
use mo_biomod     , only: k0100,k0500,k1000,k2000,k4000
use mo_bgcmean

implicit none

#include "common_clndr.h90"
#include "common_blocks.h90"

integer iogrp

integer i,j,k,l,nt
integer nhour,ny,nm,nd,dayfrac,irec(nbgcmax),cmpflg
character(len= 2) seqstring
character(len=80) fname(nbgcmax)
character(len=20) startdate
character(len=30) timeunits
real datenum,rnacc
logical append2file(nbgcmax)
data append2file /nbgcmax*.false./
save fname,irec,append2file

! set time information
timeunits=' '
startdate=' '
fname=' '
write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)')                           &
   &            'days since ',min(1800,nyear0),'-',1,'-',1,' 00:00'
write(startdate,'(i4.4,a1,i2.2,a1,i2.2,a6)')                               &
   &            nyear0,'-',nmonth0,'-',nday0,' 00:00'
datenum=time-time0-0.5*diagfq_bgc(iogrp)/nstep_in_day

! get file name
write (seqstring,'(I0.2)') nburst
if ( nmonth == 0 ) then
   call diafnm(runid,runid_len,expcnf,"hsedy."//seqstring,nstep,           &
      &         filefq_bgc(iogrp)/real(nstep_in_day),.false.,.true.,       &
      &         fname(iogrp)) ! mod_dia.F
else
   call diafnm(runid,runid_len,expcnf,"hsedm."//seqstring,nstep,           &
      &        filefq_bgc(iogrp)/real(nstep_in_day),.true.,.false.,        &
      &        fname(iogrp)) ! mod_dia.F
endif
irec(iogrp)=1
if (((fileann_bgc(iogrp) .and. nday_of_year==1 .or.                        &
   &  filemon_bgc(iogrp) .and. nday==1) .and. mod(nstep,nstep_in_day)<1)   &
   &  .or. .not.(fileann_bgc(iogrp) .or. filemon_bgc(iogrp)) .and.         &
   &  mod(nstep+.5,filefq_bgc(iogrp))<2.) then
   append2file(iogrp)=.false.
endif

! prepare output fields
if (mnproc==1) then
   write (lp,'(a,i6,a)') ' ncwrt_bgc: fields averaged over ',              &
      &                    nacc_bgc(iogrp),' steps'
   write(lp,*) 'irec(iogrp)',irec(iogrp)
endif
rnacc=1./real(nacc_bgc(iogrp))
cmpflg=GLB_COMPFLAG(iogrp)

! create output file
if (GLB_NCFORMAT(iogrp)==1) then
   call ncfopn(path1(1:path1_len)//fname(iogrp),'w','6', irec(iogrp),iotype)
elseif (GLB_NCFORMAT(iogrp)==2) then
   call ncfopn(path1(1:path1_len)//fname(iogrp),'w','h', irec(iogrp),iotype)
else
   call ncfopn(path1(1:path1_len)//fname(iogrp),'w','c', irec(iogrp),iotype)
endif

! define spatial and time dimensions
if (cmpflg/=0) then
   call ncdimc('pcomp',ip,0)
else
   call ncdims('x',itdm)
   call ncdims('y',jtdm)
endif
call ncdims('sigma',kdm)
call ncdims('depth',ddm)
call ncdims('ks',ks)
call ncdims('bounds',2)
call ncdims('time',0)
call hamoccvardef(iogrp,timeunits,calendar,cmpflg)
call nctime(datenum,calendar,timeunits,startdate)

! write auxillary dimension information
call ncwrt1('sigma','sigma',sigmar1)
call ncwrt1('depth','depth',depthslev)
call ncwrt1('depth_bnds','bounds depth',depthslev_bnds)

! Store sediment fields
call wrtsdm(jpowaic(iogrp), SDM_POWAIC(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'powdic','Porewater DIC',       ' ', 'mol C m-3')
call wrtsdm(jpowaal(iogrp), SDM_POWAAL(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'powalk','Porewater alkalinity',' ', 'eq m-3')
call wrtsdm(jpowaph(iogrp), SDM_POWAPH(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'powpho','Porewater phosphorus',' ', 'mol P m-3')
call wrtsdm(jpowaox(iogrp), SDM_POWAOX(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'powox', 'Porewater oxygen',    ' ', 'mol O2 m-3')
call wrtsdm(jpown2(iogrp),  SDM_POWN2(iogrp),  rnacc*1e3,0., cmpflg,       &
   &        'pown2', 'Porewater N2',        ' ', 'mol N2 m-3')
call wrtsdm(jpowno3(iogrp), SDM_POWNO3(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'powno3','Porewater nitrate',   ' ', 'mol N m-3')
call wrtsdm(jpowasi(iogrp), SDM_POWASI(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'powsi', 'Porewater silicate',  ' ', 'mol Si m-3')
call wrtsdm(jssso12(iogrp), SDM_SSSO12(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'ssso12','Sediment detritus',   ' ', 'mol P m-3')
call wrtsdm(jssssil(iogrp), SDM_SSSSIL(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'ssssil','Sediment silicate',   ' ', 'mol Si m-3')
call wrtsdm(jsssc12(iogrp), SDM_SSSC12(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'sssc12','Sediment CaCO3',      ' ', 'mol C m-3')
call wrtsdm(jssster(iogrp), SDM_SSSTER(iogrp), rnacc*1e3,0., cmpflg,       &
   &        'ssster','Sediment clay',       ' ', 'mol m-3')

! Store sediment burial fields
call wrtbur(jburssso12(iogrp),BUR_SSSO12(iogrp),rnacc*1e3,0.,cmpflg,       &
   & 'buro12','Burial org carbon',' ','mol P m-2')
call wrtbur(jbursssc12(iogrp),BUR_SSSC12(iogrp),rnacc*1e3,0.,cmpflg,       &
   & 'burc12','Burial calcium ',' ','mol C m-2')
call wrtbur(jburssssil(iogrp),BUR_SSSSIL(iogrp),rnacc*1e3,0.,cmpflg,       &
   & 'bursil','Burial silicate',' ','mol Si m-2')
call wrtbur(jburssster(iogrp),BUR_SSSTER(iogrp),rnacc*1e3,0.,cmpflg,       &
   & 'burter','Burial clay',' ','mol  m-2')

! close netcdf file
call ncfcls

! Initialise fields
call inisrf(jkwco2(iogrp),0.)
call inisrf(jpco2(iogrp),0.)
call inisrf(jdmsflux(iogrp),0.)
call inisrf(jco2fxd(iogrp),0.)
call inisrf(jco2fxu(iogrp),0.)
call inisrf(joxflux(iogrp),0.)
call inisrf(jniflux(iogrp),0.)
call inisrf(jn2ofx(iogrp),0.)
call inisrf(jdms(iogrp),0.)
call inisrf(jdmsprod(iogrp),0.)
call inisrf(jdms_bac(iogrp),0.)
call inisrf(jdms_uv(iogrp),0.)
call inisrf(jexport(iogrp),0.)
call inisrf(jexposi(iogrp),0.)
call inisrf(jexpoca(iogrp),0.)
call inisrf(jsrfdic(iogrp),0.)
call inisrf(jsrfalkali(iogrp),0.)
call inisrf(jsrfphosph(iogrp),0.)
call inisrf(jsrfoxygen(iogrp),0.)
call inisrf(jsrfano3(iogrp),0.)
call inisrf(jsrfsilica(iogrp),0.)
call inisrf(jsrfiron(iogrp),0.)
call inisrf(jintphosy(iogrp),0.)
call inisrf(jintnfix(iogrp),0.)
call inisrf(jintdnit(iogrp),0.)
call inisrf(jcarflx0100(iogrp),0.)
call inisrf(jcarflx0500(iogrp),0.)
call inisrf(jcarflx1000(iogrp),0.)
call inisrf(jcarflx2000(iogrp),0.)
call inisrf(jcarflx4000(iogrp),0.)
call inisrf(jcarflx_bot(iogrp),0.)
call inisrf(jbsiflx0100(iogrp),0.)
call inisrf(jbsiflx0500(iogrp),0.)
call inisrf(jbsiflx1000(iogrp),0.)
call inisrf(jbsiflx2000(iogrp),0.)
call inisrf(jbsiflx4000(iogrp),0.)
call inisrf(jbsiflx_bot(iogrp),0.)
call inisrf(jcalflx0100(iogrp),0.)
call inisrf(jcalflx0500(iogrp),0.)
call inisrf(jcalflx1000(iogrp),0.)
call inisrf(jcalflx2000(iogrp),0.)
call inisrf(jcalflx4000(iogrp),0.)
call inisrf(jcalflx_bot(iogrp),0.)
#ifdef CFC
call inisrf(jcfc11fx(iogrp),0.)
call inisrf(jcfc12fx(iogrp),0.)
call inisrf(jsf6fx(iogrp),0.)
#endif
call inisrf(jatmco2(iogrp),0.)
#ifdef DIFFAT
call inisrf(jatmo2(iogrp),0.)
call inisrf(jatmn2(iogrp),0.)
#endif
#ifdef natDIC
call inisrf(jnatco2fx(iogrp),0.)
#endif

call inilyr(jdp(iogrp),0.)
call inilyr(jdic(iogrp),0.)
call inilyr(jalkali(iogrp),0.)
call inilyr(jphosy(iogrp),0.)
call inilyr(jphosph(iogrp),0.)
call inilyr(joxygen(iogrp),0.)
call inilyr(jano3(iogrp),0.)
call inilyr(jsilica(iogrp),0.)
call inilyr(jdoc(iogrp),0.)
call inilyr(jphyto(iogrp),0.)
call inilyr(jgrazer(iogrp),0.)
call inilyr(jpoc(iogrp),0.)
call inilyr(jcalc(iogrp),0.)
call inilyr(jopal(iogrp),0.)
call inilyr(jiron(iogrp),0.)
call inilyr(jco3(iogrp),0.)
call inilyr(jph(iogrp),0.)
call inilyr(jomegac(iogrp),0.)
call inilyr(jn2o(iogrp),0.)
call inilyr(jaou(iogrp),0.)
call inilyr(jprefo2(iogrp),0.)
call inilyr(jprefpo4(iogrp),0.)
call inilyr(jprefalk(iogrp),0.)
#ifdef __c_isotopes
call inilyr(jdic13(iogrp),0.)
call inilyr(jdic14(iogrp),0.)
#endif
#ifdef AGG
call inilyr(jnos(iogrp),0.)
call inilyr(jwphy(iogrp),0.)
call inilyr(jwnos(iogrp),0.)
call inilyr(jeps(iogrp),0.)
call inilyr(jasize(iogrp),0.)
#endif
#ifdef CFC
call inilyr(jcfc11(iogrp),0.)
call inilyr(jcfc12(iogrp),0.)
call inilyr(jsf6(iogrp),0.)
#endif
#ifdef natDIC
call inilyr(jnatco3(iogrp),0.)
call inilyr(jnatalkali(iogrp),0.)
call inilyr(jnatdic(iogrp),0.)
call inilyr(jnatcalc(iogrp),0.)
call inilyr(jnatomegac(iogrp),0.)
#endif

call inilvl(jlvldic(iogrp),0.)
call inilvl(jlvlalkali(iogrp),0.)
call inilvl(jlvlphosy(iogrp),0.)
call inilvl(jlvlphosph(iogrp),0.)
call inilvl(jlvloxygen(iogrp),0.)
call inilvl(jlvlano3(iogrp),0.)
call inilvl(jlvlsilica(iogrp),0.)
call inilvl(jlvldoc(iogrp),0.)
call inilvl(jlvlphyto(iogrp),0.)
call inilvl(jlvlgrazer(iogrp),0.)
call inilvl(jlvlpoc(iogrp),0.)
call inilvl(jlvlcalc(iogrp),0.)
call inilvl(jlvlopal(iogrp),0.)
call inilvl(jlvliron(iogrp),0.)
call inilvl(jlvlco3(iogrp),0.)
call inilvl(jlvlph(iogrp),0.)
call inilvl(jlvlomegac(iogrp),0.)
call inilvl(jlvln2o(iogrp),0.)
call inilvl(jlvlaou(iogrp),0.)
call inilvl(jlvlprefo2(iogrp),0.)
call inilvl(jlvlprefpo4(iogrp),0.)
call inilvl(jlvlprefalk(iogrp),0.)
#ifdef __c_isotopes
call inilvl(jlvldic13(iogrp),0.)
call inilvl(jlvldic14(iogrp),0.)
#endif
#ifdef AGG
call inilvl(jlvlnos(iogrp),0.)
call inilvl(jlvlwphy(iogrp),0.)
call inilvl(jlvlwnos(iogrp),0.)
call inilvl(jlvleps(iogrp),0.)
call inilvl(jlvlasize(iogrp),0.)
#endif
#ifdef CFC
call inilvl(jlvlcfc11(iogrp),0.)
call inilvl(jlvlcfc12(iogrp),0.)
call inilvl(jlvlsf6(iogrp),0.)
#endif
#ifdef natDIC
call inilvl(jlvlnatco3(iogrp),0.)
call inilvl(jlvlnatalkali(iogrp),0.)
call inilvl(jlvlnatdic(iogrp),0.)
call inilvl(jlvlnatcalc(iogrp),0.)
call inilvl(jlvlnatomegac(iogrp),0.)
#endif

call inisdm(jpowaic(iogrp),0.)
call inisdm(jpowaal(iogrp),0.)
call inisdm(jpowaph(iogrp),0.)
call inisdm(jpowaox(iogrp),0.)
call inisdm(jpown2(iogrp),0.)
call inisdm(jpowno3(iogrp),0.)
call inisdm(jpowasi(iogrp),0.)
call inisdm(jssso12(iogrp),0.)
call inisdm(jssssil(iogrp),0.)
call inisdm(jsssc12(iogrp),0.)
call inisdm(jssster(iogrp),0.)

call inibur(jburssso12(iogrp),0.)
call inibur(jbursssc12(iogrp),0.)
call inibur(jburssssil(iogrp),0.)
call inibur(jburssster(iogrp),0.)

nacc_bgc(iogrp)=0

end subroutine ncwrt_onlysed

