subroutine dielectfcn(pol,epsfcn,Depsfcn)

! determines the dielectric function using an average mixing rule

use mcharge
use globals, only : ntot

implicit none
integer iR
real*8 dielPr
real*8 pol(ntot)

real*8 epsfcn(0:ntot+1)
real*8 Depsfcn(0:ntot+1)
real*8 dielW

dielW = 78.54
dielPr = dielP/dielW

do iR = 1, ntot
epsfcn(iR) = pol(iR)*dielPr + (1.0-pol(iR)) 
Depsfcn(iR) = dielPr - 1.0
enddo

epsfcn(0) = epsfcn(1) 
epsfcn(ntot+1) = epsfcn(ntot)

end subroutine
