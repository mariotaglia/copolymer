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

real*8 solvetime1, solvetime2, solveduration
integer is,ic
integer *4 ier ! Kinsol error flag
integer av1(ntot), av2(ntot)
real*8 avtmp
real*8 zc(dimR)           ! z-coordinate layer 

REAL*8 sumrhoz, meanz     ! Espesor medio pesado
real*8 pro                ! probability distribution function 
real*8 trash, trash2

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 fnorm              ! L2 norm of residual vector function fcn

external fcnelect         ! function containing the SCMFT eqs for solver
integer i,j,k,m,ii,flag,c, jj, iR,iZ ! dummy indice0s

INTEGER temp_R, temp_Z
real*8 tempr_R, tempr_Z
real*8 tmp

real*8 min1               ! variable to determine minimal position of chain
integer qqq,www,eee

integer il,inda,ncha

real*8 algo, algo2                  

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

integer in1tmp(long,2)


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

do i = 1, dimR
zc(iR)= (iR-0.5) * deltaR
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Files for multiple runs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(unit=533,file='lnq.dat')
open(unit=534,file='rgyration.dat')

end
