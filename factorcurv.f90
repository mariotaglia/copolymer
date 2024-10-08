double precision function factorcurv (bas, pos)

use globals
implicit none
integer bas, pos


select case (curvature)
case (0)
factorcurv = 1.0
case (1)
factorcurv = (float(bas+dimRini)-0.5)/(float(pos+dimRini)-0.5)
case (2)
factorcurv = ((float(bas+dimRini)-0.5)/(float(pos+dimRini)-0.5))**2
end select

return
end

double precision function jacobian(i)
use globals
use pis
use layer
implicit none
integer i
select case (curvature)
case (0)
jacobian = 1.0
case(1)
jacobian = 2.0*pi*(dfloat(i+dimRini)-0.5)*deltaR ! 2*pi*r
case(2) 
jacobian = 4.0*pi*((dfloat(i+dimRini)-0.5)**2)*deltaR**2 ! 4*pi*r^2
endselect
end 
