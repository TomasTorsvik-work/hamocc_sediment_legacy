module mo_common_bgc

use mod_xc, only: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp

! from obsolescent common_bgc.h/bgc_hamocc
real, dimension(idm,jdm)      :: bgc_dx, bgc_dy
real, dimension(idm,jdm,kdm)  :: bgc_dp, bgc_dpio
real, dimension(idm,jdm,kdm+1):: bgc_pu,bgc_pw
real, dimension(idm,jdm,kdm)  :: bgc_rho,bgc_t,bgc_s
real, dimension(idm,jdm)      :: omask
real, dimension(idm,jdm)      :: bgc_swr,bgc_fice,bgc_awnd,bgc_slp
real, dimension(idm,jdm)      :: bgc_atmco2,bgc_flxco2,bgc_flxdms

! from obsolescent common_bgc.h/bgc_hamocc_b
integer :: ldtday, ldtmonth, kpndtrun

! from obsolescent common_bgc.h/bgcc
real               :: bgcdt
integer, parameter :: nphys = 2

! from obsolescent common_bgc.h/bgc_hamocc2
real, dimension(idm,jdm) :: pglon, pglat
real                     :: bgc3dwrt, bgc2dwrt

end module mo_common_bgc
