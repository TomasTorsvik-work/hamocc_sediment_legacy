#if defined(SED_OFFLINE)
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
!use mo_bgcmean   ! FIXME: otherwise collision of GLB_FNAMETAG, ...?
use mo_sedmnt_offline ! TODO: maybe limit scope (only: GLB_FNAMETAG, ...)?
                      !       or add this routine to mo_sedmnt_offline.F90
                      !  or move namelist block from mo_sedmnt_offline.F90 to here

implicit none

#include "common_clndr.h90"
#include "common_blocks.h90"

integer, intent(in) :: iogrp

integer i,j,k,l,nt
integer nhour,ny,nm,nd,dayfrac,irec(nsedmax),cmpflg
character(len= 2) seqstring
character(len=80) fname(nsedmax)
character(len=20) startdate
character(len=30) timeunits
real datenum,rnacc
logical append2file(nsedmax)
data append2file /nsedmax*.false./
save fname,irec,append2file

! set time information
timeunits=' '
startdate=' '
fname=' '
write(timeunits,'(a11,i4.4,a1,i2.2,a1,i2.2,a6)')                           &
   &            'days since ',min(1800,nyear0),'-',1,'-',1,' 00:00'
write(startdate,'(i4.4,a1,i2.2,a1,i2.2,a6)')                               &
   &            nyear0,'-',nmonth0,'-',nday0,' 00:00'
datenum=time-time0-0.5*diagfq_sed(iogrp)/nstep_in_day

! get file name ! FIXME: we assume .not.append2file(iogrp), see ncout_hamocc.F to extend.
write (seqstring,'(I0.2)') nburst
call diafnm(runid,runid_len,expcnf,trim(GLB_FNAMETAG(iogrp))//"."//seqstring,nstep, &
   &         filefq_sed(iogrp)/real(nstep_in_day),filemon_sed(iogrp),       &
   &         fileann_sed(iogrp),fname(iogrp)) ! mod_dia.F
irec(iogrp)=1
if ( (fileann_sed(iogrp) .or. filemon_sed(iogrp))                       &
   & .or. .not.(fileann_sed(iogrp) .or. filemon_sed(iogrp)) .and.       &
   &  mod(nstep+.5,filefq_sed(iogrp))<2.) then
   append2file(iogrp) = .false.
endif   ! /FIXME

! prepare output fields
if (mnproc==1) then
   write (lp,'(a,i6,a)') ' ncwrt_sed: fields averaged over ',              &
      &                    nacc_sed(iogrp),' steps'
   write(lp,*) 'irec(iogrp)',irec(iogrp)
endif
rnacc=1./real(nacc_sed(iogrp))
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

nacc_sed(iogrp)=0

end subroutine ncwrt_onlysed
#endif
