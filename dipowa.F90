subroutine dipowa(kpie,kpje,kpke,pdlxp,pdlyp,omask,                  &
   &              bolay_)

!-----------------------------------------------------------------------
!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/dipowa.f90,v $\\
!$Revision: 1.2.20.1.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!
!-----------------------------------------------------------------------
!
! Vertical diffusion of sediment pore water tracers
!
! Copyright (C) 2001 MPI-Met, HH / Ernst Maier-Reimer
! Copyright (C) 2001 MPI-MaD, HH / S. Legutke
! Copyright (C) 2018 UiB         / Marco van Hulten
!
! Purpose:
!
! Calculate vertical diffusion of sediment pore water properties
! and diffusive flux through the ocean/sediment interface.
! Integration.
!
! Method:
!
! Implicit formulation;
! constant diffusion coefficient : 1.e-9 set in BODENSED.
! Diffusion coefficient : zcoefsu/zcoeflo for upper/lower
! sediment layer boundary.
!
! Interface:
!
!     *CALL*       *DIPOWA*
!
! Externals:
!
!     none.
!
! History:
!
! 2001-04-10   Ernst Maier-Reimer
! - Probably initial version.
! 2001-04-10   S. Legutke
! - All npowtra-1 properties are diffused in 1 go.
! JS: not mass conserving check c13/powtra/ocetra
! 2006-04-03
! - Pushed into CVS repository.
! 2018-07-13   Marco van Hulten
! - Conversion to free format etc.
! - Generalised for normal use or within sediment spin-up.
!
!-----------------------------------------------------------------------

use mo_carbch
use mo_sedmnt
use mo_biomod
use mo_param1_bgc
use mo_control_bgc

implicit none

integer, intent(in)                    :: kpie, kpje, kpke
real, dimension(kpie,kpje), intent(in) :: pdlxp, pdlyp, omask
real, dimension(kpie,kpje), intent(in) :: bolay_

integer :: i,j,k,l,iv
integer :: iv_oc                    ! index of ocetra in powtra loop

real :: sedb1(kpie,0:ks,npowtra)    ! ????
real :: zcoefsu(0:ks),zcoeflo(0:ks) ! diffusion coefficients (upper/lower)
real :: tredsy(kpie,0:kpke,3)       ! redsy for 'reduced system'?
real :: aprior                      ! start value of oceanic tracer in bottom layer
real :: bolven(kpie)                ! bottom layer ventilation rate

zcoefsu(0) = 0.0
do k = 1,ks
   ! diffusion coefficient * 1/dz     * fraction of pore water at half depths
   zcoefsu(k  ) = -sedict * seddzi(k) * porwah(k) * rdtsed
   zcoeflo(k-1) = -sedict * seddzi(k) * porwah(k) * rdtsed   ! why the same ?
enddo
zcoeflo(ks) = 0.0             ! diffusion coefficient for bottom sediment layer

!$OMP PARALLEL DO                            &
!$OMP&PRIVATE(bolven,tredsy,sedb1,aprior,iv_oc)
DO 11000 j=1,kpje

! calculate bottom ventilation rate for scaling of sediment-water exchange
!
do i = 1,kpie
   bolven(i) = 1.
enddo

k = 0
do i = 1,kpie
   tredsy(i,k,1) = zcoefsu(k)
   tredsy(i,k,3) = zcoeflo(k)
   tredsy(i,k,2) = bolven(i)*bolay_(i,j) - tredsy(i,k,1) - tredsy(i,k,3)
   !                             dz(kbo) - diff upper    - diff lower
enddo

k = 0
do iv = 1,npowtra      ! loop over pore water tracers
   iv_oc = iv
#ifdef __c_isotopes
   if (iv == ipowc13) iv_oc = isco213
   if (iv == ipowc14) iv_oc = isco214
#endif
   do i = 1,kpie
      sedb1(i,k,iv) = 0.
      if (omask(i,j) > 0.5)                                          &
         & sedb1(i,k,iv) = ocetra(i,j,kbo(i,j),iv_oc)*bolay_(i,j)*bolven(i)
   enddo
enddo

do k = 1,ks
   do i = 1,kpie
      tredsy(i,k,1) = zcoefsu(k)
      tredsy(i,k,3) = zcoeflo(k)
      tredsy(i,k,2) = seddw(k)*porwat(k) -tredsy(i,k,1) -tredsy(i,k,3)
   enddo
enddo

do iv = 1,npowtra
   do k = 1,ks
      do i = 1,kpie
         ! tracer_concentration(k[1:ks]) * porewater fraction(k) * dz(k)
         sedb1(i,k,iv) = powtra(i,j,k,iv) * porwat(k) * seddw(k)
      enddo
   enddo
enddo

do k = 1,ks
   do i = 1,kpie
      if (omask(i,j) > 0.5) then
         tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2) ! overwrites tredsy(k=0) for k=1
         !                 diff upper    / conc (k-1)
         tredsy(i,k,2)   = tredsy(i,k,2)                             &
            &             -tredsy(i,k-1,3)*tredsy(i,k,1) / tredsy(i,k-1,2)
         ! concentration  -diff lower     * diff upper   / conc(k-1)
      endif
   enddo
enddo

! diffusion from above
!
do iv = 1,npowtra
   do k = 1,ks
      do i = 1,kpie
         sedb1(i,k,iv) = sedb1(i,k,iv)                        &
            &          - tredsy(i,k-1,1) * sedb1(i,k-1,iv)
      enddo
   enddo
enddo

! sediment bottom layer
!
k = ks
do iv = 1,npowtra
   do i = 1,kpie
      if (omask(i,j) > 0.5) then
         powtra(i,j,k,iv) = sedb1(i,k,iv) / tredsy(i,k,2)
      endif
   enddo
enddo

! sediment column
!
do iv = 1,npowtra
   do k = 1,ks-1
      l = ks-k
      do i=1,kpie
          if (omask(i,j) > 0.5) then
            powtra(i,j,l,iv) = ( sedb1(i,l,iv)            &
               &             - tredsy(i,l,3) * powtra(i,j,l+1,iv) )    &
               &             / tredsy(i,l,2)
          endif
      enddo
   enddo
enddo

if (.not. lspinup_sediment) then
   !
   ! sediment ocean interface
   !
   ! NOTE: ocetra(:,:,kbo,:) cannot be replaced by ocetra_, because the former
   !       will not be updated by the latter (MvH)!
   !
   ! CAUTION - the following assumes same indices for ocetra and powtra
   do iv = 1,npowtra        ! check mo_param1_bgc.f90 for consistency
      iv_oc=iv
#ifdef __c_isotopes
      if (iv == ipowc13) iv_oc = isco213
      if (iv == ipowc14) iv_oc = isco214
#endif
      do i = 1,kpie
         l = 0
         if (omask(i,j) > 0.5) then

            aprior = ocetra(i,j,kbo(i,j),iv_oc)
            ocetra(i,j,kbo(i,j),iv_oc) =                               &
               &         ( sedb1(i,l,iv) - tredsy(i,l,3) * powtra(i,j,l+1,iv) ) &
               &         / tredsy(i,l,2)

            ! used in inventory_bgc/maschk (diagnostics)
            !
            sedfluxo(i,j,iv) = sedfluxo(i,j,iv)                          &
               &             + ocetra(i,j,kbo(i,j),iv) - aprior
#ifdef natDIC
            if (iv==isco212) ocetra(i,j,kbo(i,j),inatsco212) =                 &
               & ocetra(i,j,kbo(i,j),inatsco212) + ocetra(i,j,kbo(i,j),iv)-aprior
            if (iv==ialkali) ocetra(i,j,kbo(i,j),inatalkali) =                 &
               & ocetra(i,j,kbo(i,j),inatalkali) + ocetra(i,j,kbo(i,j),iv)-aprior
#endif
         endif
      enddo
   enddo
endif

11000 CONTINUE    ! j loop

return
end subroutine dipowa

