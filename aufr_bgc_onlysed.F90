#if defined(SED_OFFLINE)
subroutine aufr_bgc_onlysed(kpie,kpje,kpke,kplmon,rstfnm)

!-----------------------------------------------------------------------
!
! *aufr_bgc_onlysed* - read climatology to force sediment model with
!
! Based on the AUFW_BGC routine by E. Maier-Reimer (2001), S. Legutke
! (2001) and J. Schwinger (2013, 2014).  Adapted to use for reading
! bottom water properties by Marco van Hulten (2018).
!
!
!     Method:
!     -------
!     Variables for a bottom water climatology (as calculated from HAMOCC) is
!     read from a netCDF file.
!
!
!**   Interface to ocean model:
!     -------------------------
!
!     *INTEGER* *kpie*       - 1st dimension of model grid.
!     *INTEGER* *kpje*       - 2nd dimension of model grid.
!     *INTEGER* *kpke*       - 3rd (vertical) dimension of model grid.
!     *INTEGER* *kplmon*     - month of the climatology \in {1, 2, ..., 12}
!     *CHAR*    *rstfnm*     - climatology file name string
!
!-----------------------------------------------------------------------

   use netcdf
   use mo_sedmnt_offline
   use mo_control_bgc, only: io_stdo_bgc, ldtbgc, rmasko
   use mo_param1_bgc 
   use mod_xc,         only: itdm,jtdm,mnproc,xchalt
   use mod_dia,        only: iotype

   implicit none

! Function arguments
!
   integer,          intent(in)  :: kpie,kpje,kpke
   integer,          intent(in)  :: kplmon
   character(len=*), intent(in)  :: rstfnm

! Local variables
!
   integer           :: ncid,ncvarid,ncstat,ncoldmod,ncdimst(4)      &
  &                    ,nclatid,nclonid,otraid,carkid,timeid         &
  &                    ,nstart2(2),ncount2(2),nstride2(2)
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

   rmissing = rmasko

!
! Open netCDF data file
!
   IF(mnproc==1 .AND. IOTYPE==0) THEN

      write(io_stdo_bgc,*) 'Read month ', kplmon, ' from ', rstfnm
      ncstat = nf90_open(rstfnm, nf90_nowrite, ncid)
      IF ( ncstat /= NF90_NOERR ) THEN
        call xchalt('aufr_bgc_onlysed(): Problem with opening climatology')
               stop 'aufr_bgc_onlysed(): Problem with opening climatology'
      ENDIF
   ELSE IF (IOTYPE==1) THEN
#ifdef PNETCDF
      testio=1

      write(io_stdo_bgc,*) 'Read bottom seawater climatology from ', rstfnm
      write(stripestr,('(i3)')) 16
      write(stripestr2,('(i9)')) 1024*1024
      call mpi_info_create(info,ierr)
      call mpi_info_set(info,'romio_ds_read','disable',ierr)
      call mpi_info_set(info,'romio_ds_write','disable',ierr)
      call mpi_info_set(info,"striping_factor",stripestr,ierr)
      call mpi_info_set(info,"striping_unit",stripestr2,ierr)
      ncstat = NFMPI_OPEN(mpicomm,rstfnm,NF_NOWRITE,   &
                   &    info, ncid)
      IF ( ncstat /= NF_NOERR ) THEN
        call xchalt('aufr_bgc_onlysed(): Problem with opening climatology')
               stop 'aufr_bgc_onlysed(): Problem with opening climatology'
      ENDIF
#endif

      IF ( ncstat /= NF90_NOERR ) THEN
        call xchalt('(AUFW: Problem with netCDF2)')
               stop '(AUFW: Problem with netCDF2)'
      ENDIF
   ENDIF

!
! Read climatology data : ocean aquateous tracers
!--------------------------------------------------------------------
!
   CALL read_netcdf_var(ncid,'sco212',ocetra_kbo_avg(1,1,isco212),1,kplmon,iotype)
#ifdef __c_isotopes
   CALL read_netcdf_var(ncid,'sco213',ocetra_kbo_avg(1,1,isco213),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'sco214',ocetra_kbo_avg(1,1,isco214),1,kplmon,iotype)
#endif
   CALL read_netcdf_var(ncid,'alkali',ocetra_kbo_avg(1,1,ialkali),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'phosph',ocetra_kbo_avg(1,1,iphosph),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'oxygen',ocetra_kbo_avg(1,1,ioxygen),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'gasnit',ocetra_kbo_avg(1,1,igasnit),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'ano3',ocetra_kbo_avg(1,1,iano3),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'silica',ocetra_kbo_avg(1,1,isilica),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'doc',ocetra_kbo_avg(1,1,idoc),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'poc',ocetra_kbo_avg(1,1,idet),1,kplmon,iotype)
#ifdef __c_isotopes
   CALL read_netcdf_var(ncid,'poc13',ocetra_kbo_avg(1,1,idet13),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'poc14',ocetra_kbo_avg(1,1,idet14),1,kplmon,iotype)
#endif
   CALL read_netcdf_var(ncid,'phyto',ocetra_kbo_avg(1,1,iphy),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'grazer',ocetra_kbo_avg(1,1,izoo),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'calciu',ocetra_kbo_avg(1,1,icalc),1,kplmon,iotype)
#ifdef __c_isotopes
   CALL read_netcdf_var(ncid,'calciu13',ocetra_kbo_avg(1,1,icalc13),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'calciu14',ocetra_kbo_avg(1,1,icalc14),1,kplmon,iotype)
#endif
   CALL read_netcdf_var(ncid,'opal',ocetra_kbo_avg(1,1,iopal),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'n2o',ocetra_kbo_avg(1,1,ian2o),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'dms',ocetra_kbo_avg(1,1,idms),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'fdust',ocetra_kbo_avg(1,1,ifdust),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'iron',ocetra_kbo_avg(1,1,iiron),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'prefo2',ocetra_kbo_avg(1,1,iprefo2),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'prefpo4',ocetra_kbo_avg(1,1,iprefpo4),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'prefalk',ocetra_kbo_avg(1,1,iprefalk),1,kplmon,iotype)
#ifdef AGG
   CALL read_netcdf_var(ncid,'snos',ocetra_kbo_avg(1,1,inos),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'adust',ocetra_kbo_avg(1,1,iadust),1,kplmon,iotype)
#endif /*AGG*/
#ifdef ANTC14
   CALL read_netcdf_var(ncid,'antc14',ocetra_kbo_avg(1,1,iantc14),1,kplmon,iotype)
#endif
#ifdef CFC
   CALL read_netcdf_var(ncid,'cfc11',ocetra_kbo_avg(1,1,icfc11),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'cfc12',ocetra_kbo_avg(1,1,icfc12),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'sf6',ocetra_kbo_avg(1,1,isf6),1,kplmon,iotype)
#endif
#ifdef natDIC
   CALL read_netcdf_var(ncid,'natsco212',ocetra_kbo_avg(1,1,inatsco212),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'natalkali',ocetra_kbo_avg(1,1,inatalkali),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'natcalciu',ocetra_kbo_avg(1,1,inatcalc),1,kplmon,iotype)
#endif
!
! Write climatology : bottom layer depth, chemical constants and sediment variables.
!
   CALL read_netcdf_var(ncid,'bolay', bolay_avg(1,1),1,kplmon,iotype)
!  CALL read_netcdf_var(ncid,'keqb',  keqb_wrt(1,1,1),-11,kplmon)
   CALL read_netcdf_var(ncid,'K1',    keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(1,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'K2',    keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(2,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'Kb',    keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(3,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'Kw',    keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(4,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'Ks1',   keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(5,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'Kf',    keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(6,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'Ksi',   keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(7,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'K1p',   keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(8,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'K2p',   keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(9,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'K3p',   keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(10,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'Kspc',  keqb_one(1,1),1,kplmon,iotype)
   keqb_avg(11,:,:) = keqb_one(:,:)
   CALL read_netcdf_var(ncid,'prorca',prorca_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'prcaca',prcaca_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'silpro',silpro_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'produs',produs_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'bgc_t_kbo',bgc_t_kbo_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'bgc_s_kbo',bgc_s_kbo_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'bgc_rho_kbo',bgc_rho_kbo_avg(1,1),1,kplmon,iotype)
   CALL read_netcdf_var(ncid,'co3_kbo',co3_kbo_avg(1,1),1,kplmon,iotype)

   IF(mnproc==1 .AND. IOTYPE==0) THEN
      IF(mnproc==1) THEN
        ncstat = NF90_CLOSE(ncid)
        IF ( ncstat /= NF90_NOERR ) THEN
          call xchalt('(AUFW: netCDF200)')
                 stop '(AUFW: netCDF200)'
        ENDIF
      ENDIF
   ELSE IF(IOTYPE==1) THEN
#ifdef PNETCDF
      ncstat = NFMPI_CLOSE(ncid)
      IF ( ncstat /= NF_NOERR ) THEN
        call xchalt('(AUFW: PnetCDF200)')
               stop '(AUFW: PnetCDF200)'
      ENDIF
#endif
   ENDIF

   RETURN

end subroutine aufr_bgc_onlysed
#endif
