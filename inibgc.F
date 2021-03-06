      subroutine inibgc
c
      use mod_xc
      use mo_param1_bgc, only: ks,nsedtra,npowtra
      use mo_common_bgc
c
      implicit none
c
#include "param_trc.h"
#include "common_bgcs.h"
#include "common_trc.h"
#include "common_blocks.h"
#include "common_geo.h"
c
      integer i,j,k,l,nt
c  set bgc timestep to nphys*dt   
      bgcdt=nphys*baclin

c initialise 2D arrays for HAMOCC interface
      do j=1,jj
      do i=1,ii
        bgc_dx(i,j)=0.0
        bgc_dy(i,j)=0.0
        bgc_slp(i,j)=0.0
        omask (i,j)=0.0
        pglon (i,j)=plon(i,j)
        pglat (i,j)=plat(i,j)
      enddo
      enddo
c
c initialise 3D arrays for HAMOCC interface
      do k=1,kk
      do j=1,jj
      do i=1,ii
        bgc_dp  (i,j,k)=0.0
        bgc_dpio(i,j,k)=0.0
        bgc_rho (i,j,k)=0.0
        bgc_t   (i,j,k)=0.0
        bgc_s   (i,j,k)=0.0
        bgc_pu  (i,j,k)=0.0
        bgc_pw  (i,j,k)=0.0
      enddo
      enddo
      enddo
c
      do j=1,jj
      do i=1,ii
        bgc_pu  (i,j,kk+1)=0.0
        bgc_pw  (i,j,kk+1)=0.0
      enddo
      enddo
c
c initialise trc array (two time levels)
      do nt=itrbgc,itrbgc+ntrbgc-1
      do k=1,2*kk
      do j=1,jj
      do i=1,ii
        trc(i,j,k,nt)=0.0
      enddo
      enddo
      enddo
      enddo
c
c initialise sediment arrays (two time level copies)
      do j=1,jj
      do i=1,ii
        burial2(i,j,1:2,:)  = 0.0
      enddo
      enddo
c
      do k=1,2*ks
      do j=1,jj
      do i=1,ii
        sedlay2(i,j,ks,:) = 0.0
        powtra2(i,j,ks,:) = 0.0
      enddo
      enddo
      enddo
c
c initialise HAMOCC land/ocean mask
      do j=1,jj
      do l=1,isp(j)
      do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
        omask(i,j)=1.
      enddo
      enddo
      enddo
c
c
      return
      end

