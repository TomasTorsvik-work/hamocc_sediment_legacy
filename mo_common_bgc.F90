module mo_common_bgc

use mod_xc, only: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp

real, dimension(idm,jdm)      :: bgc_dx, bgc_dy
real, dimension(idm,jdm,kdm)  :: bgc_dp, bgc_dpio
real, dimension(idm,jdm,kdm+1):: bgc_pu,bgc_pw
real, dimension(idm,jdm,kdm)  :: bgc_rho,bgc_t,bgc_s
real, dimension(idm,jdm)      :: omask
real, dimension(idm,jdm)      :: bgc_swr,bgc_fice,bgc_awnd,bgc_slp
real, dimension(idm,jdm)      :: bgc_atmco2,bgc_flxco2,bgc_flxdms



end module mo_common_bgc
