subroutine aufw_bgc_onlysed(kpie,kpje,kpke,kplyear,kplmon,kplday,kpldtoce  &
   &                       ,rstfnm)

!-----------------------------------------------------------------------
!
! *aufw_bgc_onlysed* - write climatology of bottom-seawater properties
!
! Based on the AUFW_BGC routine by E. Maier-Reimer (2001), S. Legutke
! (2001) and J. Schwinger (2013, 2014).  Adapted to use for writing
! bottom water properties by Marco van Hulten (2018).
!
!
!     Method:
!     -------
!     Variables for a bottom water climatology (as calculated from HAMOCC) is
!     written to a netCDF file.
!
!
!**   Interface to ocean model:
!     -------------------------
!
!     *INTEGER* *kpie*       - 1st dimension of model grid.
!     *INTEGER* *kpje*       - 2nd dimension of model grid.
!     *INTEGER* *kpke*       - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmon*     - month of the climatology \in {1, 2, ..., 12}
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *INTEGER* *kpldtoce*   - step  in ocean restart date
!     *CHAR*    *rstfnm*     - climatology file name string
!
!-----------------------------------------------------------------------

   use netcdf
   use mo_sedmnt
   use mo_control_bgc, only: io_stdo_bgc, ldtbgc, rmasko
   use mo_param1_bgc 
   use mod_xc,         only: itdm,jtdm,mnproc,xchalt
   use mod_dia,        only: iotype

   implicit none

! Function arguments
!
   integer,          intent(in)  :: kpie,kpje,kpke
   integer,          intent(in)  :: kplyear,kplmon,kplday,kpldtoce
   character(len=*), intent(in)  :: rstfnm

! Local variables
!
   integer           :: ncid,ncvarid,ncstat,ncoldmod,ncdimst(4)      &
  &                    ,nclatid,nclonid,otraid,carkid,timeid         &
  &                    ,nstart2(2),ncount2(2),nstride2(2),idate(5)
   real              :: keqb_one(kpie,kpje)
   real              :: rmissing
#ifdef PNETCDF
   integer*4 ,save   :: info=MPI_INFO_NULL
   integer           :: mpicomm,mpierr,mpireq,mpistat
   common/xcmpii/       mpicomm,mpierr,mpireq(4),                    &
  &                     mpistat(mpi_status_size,4*max(iqr,jqr))
   save  /xcmpii/
#endif
   character(len=3)  :: stripestr
   character(len=9)  :: stripestr2
   integer           :: ierr,testio

   testio=0

   idate(1) = kplyear
   idate(2) = kplmon
   idate(3) = kplday
   idate(4) = kpldtoce
   idate(5) = ldtbgc
   IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Writing climatology at date :'              &
     &           ,'YY=',idate(1),' MM=',idate(2),' day=',idate(3)
      WRITE(io_stdo_bgc,*) 'Ocean model step number is ',idate(4)
      WRITE(io_stdo_bgc,*) 'Bgc   model step number is ',idate(5)
   ENDIF

   rmissing = rmasko

!
! Open netCDF data file
!
   IF(mnproc==1 .AND. IOTYPE==0) THEN

      write(io_stdo_bgc,*) 'Write bottom seawater climatology to ', rstfnm
      if (kplmon == 1) then
         ncstat = nf90_create(rstfnm, nf90_64bit_offset, ncid)
      else
         ncstat = nf90_open(rstfnm, nf90_write, ncid)
      endif
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('aufw_bgc_onlysed(): Problem with opening climatology')
               stop 'aufw_bgc_onlysed(): Problem with opening climatology'
      ENDIF
   ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      testio=1

      write(io_stdo_bgc,*) 'Write bottom seawater climatology to ', rstfnm
      write(stripestr,('(i3)')) 16
      write(stripestr2,('(i9)')) 1024*1024
      call mpi_info_create(info,ierr)
      call mpi_info_set(info,'romio_ds_read','disable',ierr)
      call mpi_info_set(info,'romio_ds_write','disable',ierr)
      call mpi_info_set(info,"striping_factor",stripestr,ierr)
      call mpi_info_set(info,"striping_unit",stripestr2,ierr)
      ncstat = NFMPI_CREATE(mpicomm,rstfnm,             &
     &       IOR(nf_clobber,nf_64bit_offset),info,ncid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('aufw_bgc_onlysed(): Problem with opening climatology')
               stop 'aufw_bgc_onlysed(): Problem with opening climatology'
      ENDIF
#endif

      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF2)')
               stop '(AUFW: Problem with netCDF2)'
      ENDIF
   ENDIF

!
! Define dimension
! ----------------------------------------------------------------------    
!
if (kplmon == 1) then ! only for first month
   IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_DEF_DIM(ncid, 'lon', itdm, nclonid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF2)')
               stop '(AUFW: Problem with netCDF2)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'lat', jtdm, nclatid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF3)')
               stop '(AUFW: Problem with netCDF3)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'otra', nocetra, otraid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF7)')
               stop '(AUFW: Problem with netCDF7)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'cark', 11, carkid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF7)')
               stop '(AUFW: Problem with netCDF7)'
      ENDIF

      ncstat = NF90_DEF_DIM(ncid, 'time', 12, timeid)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF8)')
               stop '(AUFW: Problem with netCDF8)'
      ENDIF

   ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_DEF_DIM(ncid, 'lon', itdm, nclonid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF2)')
               stop '(AUFW: Problem with PnetCDF2)'
      ENDIF

      ncstat = NFMPI_DEF_DIM(ncid, 'lat', jtdm, nclatid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF3)')
               stop '(AUFW: Problem with PnetCDF3)'
      ENDIF

      ncstat = NFMPI_DEF_DIM(ncid, 'otra', nocetra, otraid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF5)')
               stop '(AUFW: Problem with PnetCDF5)'
      ENDIF

      ncstat = NFMPI_DEF_DIM(ncid, 'cark', 11, carkid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF6)')
               stop '(AUFW: Problem with PnetCDF6)'
      ENDIF

      ncstat = NFMPI_DEF_DIM(ncid, 'time', 12, timeid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF8)')
               stop '(AUFW: Problem with PnetCDF8)'
      ENDIF
#endif
   ENDIF
endif ! only for first month

!
! Define global attributes
! ----------------------------------------------------------------------    
!
if (kplmon == 1) then ! only for first month
   IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'title'               &
     &, 'Forcing data for the stand-alone sediment module') 
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'history'             &
     &, 'Forcing data for the stand-alone sediment module')
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'conventions'         &
     &,'COARDS')
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF9)')
               stop '(AUFW: Problem with netCDF9)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL,'source'              &
     &, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF10)')
               stop '(AUFW: Problem with netCDF10)'
      ENDIF

      ncstat = NF90_PUT_ATT(ncid, NF90_GLOBAL, 'date', idate)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF11)')
               stop '(AUFW: Problem with netCDF11)'
      ENDIF

   ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      clen=len('Forcing data for the stand-alone sediment module')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'title'             &
     &, clen,'Forcing data for the stand-alone sediment module') 
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF9)')
               stop '(AUFW: Problem with PnetCDF9)'
      ENDIF
      clen=len('Forcing data for the stand-alone sediment module')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'history'           &
     &, clen,'Forcing data for the stand-alone sediment module')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF9)')
               stop '(AUFW: Problem with PnetCDF9)'
      ENDIF
      clen=6
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'conventions'       &
     &,clen, 'COARDS')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF9)')
               stop '(AUFW: Problem with PnetCDF9)'
      ENDIF
      clen=len('Marine bgc model output HOPC68/grob')
      ncstat = NFMPI_PUT_ATT_TEXT(ncid, NF_GLOBAL,'source'            &
     &,clen, 'Marine bgc model output HOPC68/grob')
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF10)')
               stop '(AUFW: Problem with PnetCDF10)'
      ENDIF
      clen=5
      ncstat = NFMPI_PUT_ATT_INT(ncid, NF_GLOBAL, 'date',              &
     &                           nf_int, clen, idate)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF11)')
               stop '(AUFW: Problem with netCDF11)'

      ENDIF
#endif
   ENDIF
endif ! only for first month

! 
! Define variables : advected ocean tracer
! ----------------------------------------------------------------------    
!
if (kplmon == 1) then ! only for first month
      ncdimst(1) = nclonid
      ncdimst(2) = nclatid
      ncdimst(3) = timeid
      ncdimst(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,6,'sco212',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',13, 'Dissolved CO2',rmissing,22,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,6,'sco213',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',15, 'Dissolved CO213',rmissing,22,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'sco214',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',15, 'Dissolved CO214',rmissing,22,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,6,'alkali',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',10,'Alkalinity',rmissing,25,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'phosph',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',19,'Dissolved phosphate',rmissing,28,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'oxygen',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',16,'Dissolved oxygen',                             &
          rmissing,31,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'gasnit',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',21,'Gaseous nitrogen (N2)',                        &
          rmissing,34,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'ano3',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',17,'Dissolved nitrate',                            &
          rmissing,34,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'silica',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',22,'Silicid acid (Si(OH)4)',                       &
          rmissing,40,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'doc',3,ncdimst,ncvarid,              &
     &    6,'mol/kg',24,'Dissolved organic carbon',                     &
     &    rmissing,40,io_stdo_bgc) 

      CALL NETCDF_DEF_VARDB(ncid,3,'poc',3,ncdimst,ncvarid,              &
     &    6,'mol/kg',25,'Particulate organic carbon',                   &
     &    rmissing,46,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,5,'poc13',3,ncdimst,ncvarid,            &
     &    7,'molC/kg',28,'Particulate organic carbon13',                &
     &    rmissing,46,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'poc14',3,ncdimst,ncvarid,            &
     &    7,'molC/kg',28,'Particulate organic carbon14',                &
     &    rmissing,46,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,5,'phyto',3,ncdimst,ncvarid,            &
     &    7,'molP/kg',27,'Phytoplankton concentration',                 &
     &    rmissing,28,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'grazer',3,ncdimst,ncvarid,           &
     &    7,'molP/kg',25,'Zooplankton concentration',                   &
     &    rmissing,29,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'calciu',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',17,'Calcium carbonate',                            &
     &    rmissing,30,io_stdo_bgc)

#ifdef __c_isotopes
      CALL NETCDF_DEF_VARDB(ncid,8,'calciu13',3,ncdimst,ncvarid,         &
     &    7,'molC/kg',19,'Calcium carbonate13',                         &
     &    rmissing,30,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,8,'calciu14',3,ncdimst,ncvarid,         &
     &    7,'molC/kg',19,'Calcium carbonate14',                         &
     &    rmissing,30,io_stdo_bgc)
#endif

      CALL NETCDF_DEF_VARDB(ncid,4,'opal',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',15,'Biogenic silica',                              &
     &    rmissing,31,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'n2o',3,ncdimst,ncvarid,              &
     &    6,'mol/kg',12,'laughing gas',                                 &
     &    rmissing,32,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'dms',3,ncdimst,ncvarid,              &
     &    6,'mol/kg',15 ,'DiMethylSulfide',                             &
     &    rmissing,33,io_stdo_bgc)            

      CALL NETCDF_DEF_VARDB(ncid,5,'fdust',3,ncdimst,ncvarid,            &
     &    5,'kg/kg',19,'Non-aggregated dust',                           &
     &    rmissing,34,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'iron',3,ncdimst,ncvarid,             &
     &    6,'mol/kg',14,'Dissolved iron',                               &
     &    rmissing,35,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'prefo2',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',16,'Preformed oxygen',                             &
          rmissing,43,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'prefpo4',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',19,'Preformed phosphate',                          &
          rmissing,43,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'prefalk',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',20,'Preformed alkalinity',                         &
          rmissing,43,io_stdo_bgc)

#ifdef AGG
      CALL NETCDF_DEF_VARDB(ncid,4,'snos',3,ncdimst,ncvarid,             &
     &    3,'1/g',38,'marine snow aggregates per g sea water',          &
     &    rmissing,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'adust',3,ncdimst,ncvarid,            &
     &    4,'g/kg',15,'Aggregated dust',                                &
     &    rmissing,42,io_stdo_bgc)
#endif /*AGG*/   

#ifdef ANTC14
      CALL NETCDF_DEF_VARDB(ncid,6,'antc14',3,ncdimst,ncvarid,           &
     &    6,'mol/kg',17,'anthropogenic C14',                            &
     &    rmissing,41,io_stdo_bgc)
#endif
#ifdef CFC
      CALL NETCDF_DEF_VARDB(ncid,5,'cfc11',3,ncdimst,ncvarid,            &
     &    6,'mol/kg',5,'CFC11',                                         &
     &    rmissing,41,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,5,'cfc12',3,ncdimst,ncvarid,            &
     &    6,'mol/kg',5,'CFC12',                                         &
     &    rmissing,41,io_stdo_bgc)     

      CALL NETCDF_DEF_VARDB(ncid,3,'sf6',3,ncdimst,ncvarid,              &
     &    6,'mol/kg',4,'SF-6',                                          &
     &    rmissing,41,io_stdo_bgc)     
#endif
#ifdef natDIC
      CALL NETCDF_DEF_VARDB(ncid,9,'natsco212',3,ncdimst,ncvarid,       &
     &   6,'mol/kg',21, 'Natural dissolved CO2',rmissing,22,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,9,'natalkali',3,ncdimst,ncvarid,       &
     &    6,'mol/kg',18,'Natural alkalinity',rmissing,25,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,9,'natcalciu',3,ncdimst,ncvarid,       &
     &    6,'mol/kg',25,'Natural calcium carbonate',                    &
     &    rmissing,30,io_stdo_bgc)
#endif

! 
! Define variables : chemical equilibrium constants
! ----------------------------------------------------------------------    
!
      ncdimst(1) = nclonid
      ncdimst(2) = nclatid
      ncdimst(3) = timeid
      ncdimst(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,2,'K1',3,ncdimst,ncvarid,          &
     &    1,'M',17,'[H][HCO3]/[H2CO3]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,2,'K2',3,ncdimst,ncvarid,          &
     &    1,'M',15,'[H][CO3]/[HCO3]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,2,'Kb',3,ncdimst,ncvarid,          &
     &    1,'M',15,'[H][BO2]/[HBO2]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,2,'Kw',3,ncdimst,ncvarid,          &
     &    3,'M^2',7,'[H][OH]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'Ks1',3,ncdimst,ncvarid,          &
     &    1,'M',15,'[H][SO4]/[HSO4]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,2,'Kf',3,ncdimst,ncvarid,          &
     &    1,'M',11,'[H][F]/[HF]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'Ksi',3,ncdimst,ncvarid,          &
     &    1,'M',23,'[H][SiO(OH)3]/[Si(OH)4]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'K1p',3,ncdimst,ncvarid,          &
     &    1,'M',18,'[H][H2PO4]/[H3PO4]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'K2p',3,ncdimst,ncvarid,          &
     &    1,'M',17,'[H][HPO4]/[H2PO4]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,3,'K3p',3,ncdimst,ncvarid,          &
     &    1,'M',15,'[H][PO4]/[HPO4]',      &
     &    rmissing,69,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,4,'Kspc',3,ncdimst,ncvarid,          &
     &    3,'M^2',16,'[Ca2+]_T [CO3]_T',      &
     &    rmissing,69,io_stdo_bgc)

! 
! Define variables : bottom layer thickness; sedimentation variables
! ----------------------------------------------------------------------    
!
      ncdimst(1) = nclonid
      ncdimst(2) = nclatid
      ncdimst(3) = timeid
      ncdimst(4) = 0

      CALL NETCDF_DEF_VARDB(ncid,5,'bolay',3,ncdimst,ncvarid,          &
     &    1,'m',31,'Thickness of the bottom gridbox',       &
     &    rmissing,72,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'prorca',3,ncdimst,ncvarid,          &
     &    10,'kmol/m^2/s',23,'Sedimentation of carbon',       &
     &    rmissing,77,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'prcaca',3,ncdimst,ncvarid,          &
     &    10,'kmol/m^2/s',34,'Sedimentation of calcium carbonate',       &
     &    rmissing,77,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'silpro',3,ncdimst,ncvarid,          &
     &    10,'kmol/m^2/s',24,'Sedimentation of silicon',       &
     &    rmissing,77,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,6,'produs',3,ncdimst,ncvarid,          &
     &    10,'kmol/m^2/s',21,'Sedimentation of dust',       &
     &    rmissing,77,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,9,'bgc_t_kbo',3,ncdimst,ncvarid,          &
     &    6,'degr C',33,'Temperature in the bottom gridbox',       &
     &    rmissing,73,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,9,'bgc_s_kbo',3,ncdimst,ncvarid,          &
     &    3,'psu',30,'Salinity in the bottom gridbox',       &
     &    rmissing,74,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,11,'bgc_rho_kbo',3,ncdimst,ncvarid,          &
     &    6,'g/cm^3',29,'Density in the bottom gridbox',       &
     &    rmissing,75,io_stdo_bgc)

      CALL NETCDF_DEF_VARDB(ncid,7,'co3_kbo',3,ncdimst,ncvarid,          &
     &    6,'mol/kg',41,'Dissolved carbonate in the bottom gridbox',       &
     &    rmissing,76,io_stdo_bgc)

   IF(mnproc==1 .AND. IOTYPE==0) THEN
      ncstat = NF90_ENDDEF(ncid)
      
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF00)')
               stop '(AUFW: Problem with netCDF00)'
      ENDIF

!
! Set fill mode
! ----------------------------------------------------------------------    
!
      ncstat = NF90_SET_FILL(ncid,NF90_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF97)')
               stop '(AUFW: Problem with netCDF97)'
      ENDIF

   ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_ENDDEF(ncid)

      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: Problem with PnetCDF00)')
               stop '(AUFW: Problem with PnetCDF00)'
      ENDIF
#endif
   ENDIF
endif ! only for first month

!
! Write climatology data : ocean aquateous tracers
!--------------------------------------------------------------------
!
   CALL write_netcdf_var(ncid,'sco212',ocetra_kbo_avg(1,1,isco212),1,kplmon)
#ifdef __c_isotopes
   CALL write_netcdf_var(ncid,'sco213',ocetra_kbo_avg(1,1,isco213),1,kplmon)
   CALL write_netcdf_var(ncid,'sco214',ocetra_kbo_avg(1,1,isco214),1,kplmon)
#endif
   CALL write_netcdf_var(ncid,'alkali',ocetra_kbo_avg(1,1,ialkali),1,kplmon)
   CALL write_netcdf_var(ncid,'phosph',ocetra_kbo_avg(1,1,iphosph),1,kplmon)
   CALL write_netcdf_var(ncid,'oxygen',ocetra_kbo_avg(1,1,ioxygen),1,kplmon)
   CALL write_netcdf_var(ncid,'gasnit',ocetra_kbo_avg(1,1,igasnit),1,kplmon)
   CALL write_netcdf_var(ncid,'ano3',ocetra_kbo_avg(1,1,iano3),1,kplmon)
   CALL write_netcdf_var(ncid,'silica',ocetra_kbo_avg(1,1,isilica),1,kplmon)
   CALL write_netcdf_var(ncid,'doc',ocetra_kbo_avg(1,1,idoc),1,kplmon)
   CALL write_netcdf_var(ncid,'poc',ocetra_kbo_avg(1,1,idet),1,kplmon)
#ifdef __c_isotopes
   CALL write_netcdf_var(ncid,'poc13',ocetra_kbo_avg(1,1,idet13),1,kplmon)
   CALL write_netcdf_var(ncid,'poc14',ocetra_kbo_avg(1,1,idet14),1,kplmon)
#endif
   CALL write_netcdf_var(ncid,'phyto',ocetra_kbo_avg(1,1,iphy),1,kplmon)
   CALL write_netcdf_var(ncid,'grazer',ocetra_kbo_avg(1,1,izoo),1,kplmon)
   CALL write_netcdf_var(ncid,'calciu',ocetra_kbo_avg(1,1,icalc),1,kplmon)
#ifdef __c_isotopes
   CALL write_netcdf_var(ncid,'calciu13',ocetra_kbo_avg(1,1,icalc13),1,kplmon)
   CALL write_netcdf_var(ncid,'calciu14',ocetra_kbo_avg(1,1,icalc14),1,kplmon)
#endif
   CALL write_netcdf_var(ncid,'opal',ocetra_kbo_avg(1,1,iopal),1,kplmon)
   CALL write_netcdf_var(ncid,'n2o',ocetra_kbo_avg(1,1,ian2o),1,kplmon)
   CALL write_netcdf_var(ncid,'dms',ocetra_kbo_avg(1,1,idms),1,kplmon)
   CALL write_netcdf_var(ncid,'fdust',ocetra_kbo_avg(1,1,ifdust),1,kplmon)
   CALL write_netcdf_var(ncid,'iron',ocetra_kbo_avg(1,1,iiron),1,kplmon)
   CALL write_netcdf_var(ncid,'prefo2',ocetra_kbo_avg(1,1,iprefo2),1,kplmon)
   CALL write_netcdf_var(ncid,'prefpo4',ocetra_kbo_avg(1,1,iprefpo4),1,kplmon)
   CALL write_netcdf_var(ncid,'prefalk',ocetra_kbo_avg(1,1,iprefalk),1,kplmon)
#ifdef AGG
   CALL write_netcdf_var(ncid,'snos',ocetra_kbo_avg(1,1,inos),1,kplmon)
   CALL write_netcdf_var(ncid,'adust',ocetra_kbo_avg(1,1,iadust),1,kplmon)
#endif /*AGG*/
#ifdef ANTC14
   CALL write_netcdf_var(ncid,'antc14',ocetra_kbo_avg(1,1,iantc14),1,kplmon)
#endif
#ifdef CFC
   CALL write_netcdf_var(ncid,'cfc11',ocetra_kbo_avg(1,1,icfc11),1,kplmon)
   CALL write_netcdf_var(ncid,'cfc12',ocetra_kbo_avg(1,1,icfc12),1,kplmon)
   CALL write_netcdf_var(ncid,'sf6',ocetra_kbo_avg(1,1,isf6),1,kplmon)
#endif
#ifdef natDIC
   CALL write_netcdf_var(ncid,'natsco212',ocetra_kbo_avg(1,1,inatsco212),1,kplmon)
   CALL write_netcdf_var(ncid,'natalkali',ocetra_kbo_avg(1,1,inatalkali),1,kplmon)
   CALL write_netcdf_var(ncid,'natcalciu',ocetra_kbo_avg(1,1,inatcalc),1,kplmon)
#endif
!
! Write climatology : bottom layer depth, chemical constants and sediment variables.
!
   CALL write_netcdf_var(ncid,'bolay', bolay_avg(1,1),1,kplmon)
!  CALL write_netcdf_var(ncid,'keqb',  keqb_wrt(1,1,1),-11,kplmon)
   keqb_one = keqb_avg(1,:,:)
   CALL write_netcdf_var(ncid,'K1',    keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(2,:,:)
   CALL write_netcdf_var(ncid,'K2',    keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(3,:,:)
   CALL write_netcdf_var(ncid,'Kb',    keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(4,:,:)
   CALL write_netcdf_var(ncid,'Kw',    keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(5,:,:)
   CALL write_netcdf_var(ncid,'Ks1',   keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(6,:,:)
   CALL write_netcdf_var(ncid,'Kf',    keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(7,:,:)
   CALL write_netcdf_var(ncid,'Ksi',   keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(8,:,:)
   CALL write_netcdf_var(ncid,'K1p',   keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(9,:,:)
   CALL write_netcdf_var(ncid,'K2p',   keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(10,:,:)
   CALL write_netcdf_var(ncid,'K3p',   keqb_one(1,1),1,kplmon)
   keqb_one = keqb_avg(11,:,:)
   CALL write_netcdf_var(ncid,'Kspc',  keqb_one(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'prorca',prorca_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'prcaca',prcaca_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'silpro',silpro_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'produs',produs_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'bgc_t_kbo',bgc_t_kbo_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'bgc_s_kbo',bgc_s_kbo_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'bgc_rho_kbo',bgc_rho_kbo_avg(1,1),1,kplmon)
   CALL write_netcdf_var(ncid,'co3_kbo',co3_kbo_avg(1,1),1,kplmon)

   IF(mnproc==1 .AND. IOTYPE==0) THEN
      IF(mnproc==1) THEN
        ncstat = NF90_CLOSE(ncid)
        IF ( ncstat .NE. NF90_NOERR ) THEN
          call xchalt('(AUFW: netCDF200)')
                 stop '(AUFW: netCDF200)'
        ENDIF
      ENDIF
   ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_CLOSE(ncid)
      IF ( ncstat .NE. NF_NOERR ) THEN
        call xchalt('(AUFW: PnetCDF200)')
               stop '(AUFW: PnetCDF200)'
      ENDIF
#endif
   ENDIF
   IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*) 'End of AUFW_BGC'
      WRITE(io_stdo_bgc,*) '***************'
   ENDIF

   RETURN
end
