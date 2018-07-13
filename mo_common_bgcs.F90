module mo_common_bgcs

use mod_xc, only: idm, jdm
use mo_param1_bgc

real, dimension(idm,jdm,2*ks,nsedtra)  :: sedlay2
real, dimension(idm,jdm,2*ks,npowtra)  :: powtra2
real, dimension(idm,jdm,2,nsedtra)     :: burial2

end module mo_common_bgcs
