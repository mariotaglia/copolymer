subroutine dielectfcn(pol,epsfcn,Depsfcn)

! determines the dielectric function using an average mixing rule

use mcharge
use globals, only : dimR,dimZ

implicit none
integer iR,iZ
real*8 dielPr
real*8 pol(dimR,dimZ)

real*8 epsfcn(0:dimR+1,dimZ)
real*8 Depsfcn(0:dimR+1,dimZ)
real*8 dielW

dielW = 78.54
dielPr = dielP/dielW

do iR = 1, dimR
do iZ = 1, dimZ
epsfcn(iR,iZ) = pol(iR,iZ)*dielPr + (1.0-pol(iR,iZ)) 
Depsfcn(iR,iZ) = dielPr - 1.0
enddo
enddo

!epsfcn(0) = epsfcn(1) 
!epsfcn(ntot+1) = epsfcn(ntot) OJO boundary conditions  

end subroutine
