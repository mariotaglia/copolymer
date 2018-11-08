subroutine init

use globals
use mcharge 
use layer
use volume
use bulk
use longs
use MPI
use pis
implicit none

real*8 zc(dimR)           ! z-coordinate layer 

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 

external fcnelect         ! function containing the SCMFT eqs for solver
integer i, iR ! dummy indice0s

! MPI
integer tag
parameter(tag = 0)

error = 1e-6

!     initializations of variables 

pi=dacos(-1.0d0)          ! pi = arccos(-1) 
itmax=200                 ! maximum number of iterations       
n=ntot                    ! size of lattice

vsol=0.030                ! volume solvent molecule in (nm)^3
vpol(:)=vpol(:)/vsol  ! volume polymer segment in units of vsol
vpol_a(:)=vpol_a(:)/vsol
vpol_b(:)=vpol_b(:)/vsol

vchain=0.0
do i=1,long
  vchain=vchain+vpol(segpoorsv(i))
enddo

vneg=4/3*pi*0.2**3/vsol !volume of anion in units of vsol
vpos=4/3*pi*0.2**3/vsol !volume of cation in units of vsol 

pKw=14.0

cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
rhosalt=Csalt*Na/(1.0d24) !salt conc. in unit of nº of particles/nm³
wperm = 0.114 !water permitivity in units of e^2/kT.nm

if(pHbulk.le.7) then  ! pH<= 7
  xposbulk= rhosalt*vsol*vpos
  xnegbulk= rhosalt*vsol*vneg + (xHplusbulk-xOHminbulk)*vneg ! NaCl+ HCl  
else                  ! pH >7 
  xposbulk= rhosalt*vsol*vpos + (xOHminbulk-xHplusbulk)*vpos ! NaCl+ NaOH   
  xnegbulk= rhosalt*vsol*vneg
endif


xsolbulk=1-xposbulk-xnegbulk-xHplusbulk-xOHminbulk ! bulk volume fraction of solvent 

if (Nacids.gt.0) then
  do i=1,Nacids
    Ka(i) = 10**(-pKa(i))
    Ka(i) = (Ka(i)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant
  enddo
endif

if (Nbasics.gt.0) then
  do i=1,Nbasics
    Kb(i) = 10**(-pKb(i))
    Kb(i) = (Kb(i)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant
  enddo
endif

expmupos=xposbulk/vpos/xsolbulk**vpos  
expmuneg=xnegbulk/vneg/xsolbulk**vneg           
expmuHplus=xHplusbulk/xsolbulk ! vHplus=vsol
expmuOHmin=xOHminbulk/xsolbulk ! vOHminus=vsol

do iR = 1, dimR
zc(iR)= (iR-0.5) * deltaR
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Files for multiple runs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(unit=533,file='lnq.dat')
open(unit=534,file='rgyration.dat')

end
