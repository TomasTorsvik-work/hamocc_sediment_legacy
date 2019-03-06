subroutine read_netcdf_var(ncid,desc,arr,klev,time,typeio)

!-----------------------------------------------------------------------
!
! *read_netcdf_var* - read variable as netCDF data
!
! Based on write_netcdf_var().               TODO: description
! This will be done for all levels klev, and for the given time index.
! The first index is always idm and second jdm, so fix the first index
! of keqb(:,1:idm,1:jdm,:) when using this routine.
! The netCDF file is only accessed by mnproc=1.
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
real, intent(inout)          :: arr(idm,jdm,klev)
integer, intent(in)          :: typeio

integer :: ndims
real :: arr_g(itdm,jtdm)
real, allocatable :: arr_l(:,:)
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
     !call xcaget(arr_g,arr_l,1)
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
         ncstat = nf90_get_var(ncid,ncvarid,arr_g,start,count)
         if (ncstat /= nf90_noerr) then
            write(lp,'(4a)') 'nf90_get_var: ',trim(desc),': ', &
               &              nf90_strerror(ncstat)
            call xchalt('(write_netcdf_var)')
                   stop '(write_netcdf_var)'
         endif
      endif
      call xcaput(arr_g,arr_l,1)
      do j=1,jdm
         do i=1,idm
            arr(i,j,k)=arr_l(i,j)
         enddo
      enddo
   enddo

   deallocate(start,count)
else
   call xchalt('PNETCDF/iotype==1 removed from write_netcdf_var().')
          stop 'PNETCDF/iotype==1 removed from write_netcdf_var().'
endif

END
