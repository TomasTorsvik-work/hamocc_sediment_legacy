subroutine write_netcdf_var(ncid,desc,arr,klev,time)

!-----------------------------------------------------------------------
!
! MPI gather and write variable as netCDF data
!
! Copyright (C) 2018 Marco van Hulten <Marco.Hulten@uib.no> et al.
!                    Geophysical Institute @ University of Bergen
!
! This subroutine is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Method:
!
! Gathers a global variable from all PEs and writes it to a netCDF file.
! This will be done for all levels klev, and for the given time index.
! The first index is always idm and second jdm, so fix the first index
! of keqb(:,1:idm,1:jdm,:) when using this routine.
! The netCDF file is only accessed by mnproc=1.
!
! History:
!
! Rewritten by Marco van Hulten (2018).  Improved logic and readability,
! removed PNETCDF code.
!
!-----------------------------------------------------------------------

use netcdf
use mo_control_bgc
use mod_xc
use mod_dia, only : iotype      

implicit none

#include <mpif.h>

integer, intent(in)          :: ncid, klev, time
character(len=*), intent(in) :: desc
real, intent(in)             :: arr(idm, jdm, klev)

integer :: ndims
real :: arr_g(itdm,jtdm)
real, allocatable :: arr_g1(:,:,:), arr_l(:,:)
integer ncstat, ncvarid, k, i, j
integer, allocatable :: start(:), count(:) 

ndims = 2
if (klev > 1) then
   ndims = ndims + 1
endif
if (time > 0) then
   ndims = ndims + 1
endif

if (iotype == 0) then
   allocate(start(ndims))
   allocate(count(ndims))
   allocate(arr_l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
   arr_l=0.0
   start(1) = 1
   count(1) = itdm
   start(2) = 1
   count(2) = jtdm
   if (klev > 1) then
      start(3) = 1
      count(3) = 1
   endif
   if (time > 0) then
      start(ndims) = time
      count(ndims) = 1
   endif

   do k = 1, klev
      do j = 1, jdm
         do i = 1, idm
            arr_l(i,j) = arr(i,j,k)
         enddo
      enddo
      call xcaget(arr_g,arr_l,1)
      if (mnproc == 1) then
         if (k > 1) then
            start(3) = k
            count(3) = 1
         endif
         ncstat = nf90_inq_varid(ncid,desc,ncvarid)
         if (ncstat /= nf90_noerr) then
            write(lp,'(4a)') 'nf90_inq_varid: ',trim(desc),': ', &
               &              nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
         endif
         ncstat = nf90_put_var(ncid,ncvarid,arr_g,start,count)
         if (ncstat /= nf90_noerr) then
            write(lp,'(4a)') 'nf90_put_var: ',trim(desc),': ', &
               &              nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
         endif
      endif
   enddo

   deallocate(start,count)
else if (iotype == 1) then
#ifdef PNETCDF
   allocate(istart(ndims))
   allocate(icount(ndims))
   allocate(arr_l(ii,jj,klev))
   arr_l=0.0
   istart(1) = 1
   icount(1) = itdm
   istart(2) = 1
   icount(2) = jtdm
   if (klev > 1) then
      istart(3) = 1
      icount(3) = klev
   endif
   if (time > 0) then
      istart(ndims) = time
      icount(ndims) = 1
   endif

   istart(1)=1
   istart(2)=j0+1

   if (mproc == mpe_1(nproc)) then
      icount(1) = itdm
      icount(2) = jj
   else
      do i = 1, ndims
         icount(i)=0
      enddo
   endif

   do k = 1, klev
      do j = 1, jdm
         do i = 1, idm
            arr_l(i,j) = arr(i,j,k)
         enddo
      enddo
      allocate(arr_g1(itdm,jj,klev))
      arr_g1 = 0.0
      call xcgetrow(arr_g1, arr_l, klev)

      ncstat = nfmpi_inq_varid(ncid,desc,ncvarid)
      if (ncstat /= nf_noerr) then
         write(lp,'(4a)') 'nfmpi_inq_varid: ',trim(desc),': ', &
               &             nfmpi_strerror(ncstat)
            call xchalt('(write_pnetcdf_var)')
                   stop '(write_pnetcdf_var)'
         endif
         ncstat=nfmpi_put_vara_double_all(ncid,ncvarid,istart, &
               &                          icount,arr_g1)
         if (ncstat /= nf_noerr) then
            write(lp,'(4a)') 'nfmpi_put_var: ',trim(desc),': ', &
               &             nfmpi_strerror(ncstat)
            call xchalt('(write_pnetcdf_var)')
                   stop '(write_pnetcdf_var)'
         endif
      endif
   enddo

   deallocate(istart,icount,arr_g1)
!   call xchalt('PNETCDF/iotype==1 removed from write_netcdf_var().')
!          stop 'PNETCDF/iotype==1 removed from write_netcdf_var().'
#endif
endif

END
