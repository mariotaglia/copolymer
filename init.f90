subroutine init

use globals
use mcharge 
use layer
use volume
use bulk
use longs
use MPI
use pis
use mkai
implicit none


integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 

external fcnelect         ! function containing the SCMFT eqs for solver
integer i, iR, ia, ib, is ! dummy indices
integer NC
character*10 lnqfile, rogfile
real*8 chargebalance 

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


vneg=4/3*pi*r_neg**3/vsol !volume of anion in units of vsol
vpos=4/3*pi*r_pos**3/vsol !volume of cation in units of vsol 

pKw=14.0

wperm = 0.114 !water permitivity in units of e^2/kT.nm

!! Bulk number density and volume fraction of species !!

cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
rhosalt = Csalt*Na/(1.0d24) !salt conc. in unit of nº of particles/nm³

rhopolbulk = 0.
xpolbulk = 0.
avpolbulk = 0.

!! volume fraction of GC components and their beads

do NC = 1, Ncomp
   if (flagGC(NC).eq.1) then
     rhopolbulk(NC) = Cpolbulk(NC)*Na/(1.0d24) ! bulk conc. in units of n of particles/nm³ 
     do i =1,long(NC)
       is = segpoorsv(i,NC)
       avpolbulk(is,NC) = avpolbulk(is,NC) + rhopolbulk(NC) * vpol(is) * vsol 
       xpolbulk(NC) = xpolbulk(NC) + rhopolbulk(NC) * vpol(is) * vsol 
     enddo
   endif
enddo   


!! State of charge of titrable beads in the bulk !!

fAmin_bulk = 0.
fBHplus_bulk = 0.

if (Nacids.gt.0) then
   do i=1,Nacids
      Ka(i) = 10**(-pKa(i))
      fAmin_bulk(i) =  Ka(i)/cHplus / (1.0 + Ka(i)/cHplus)
   enddo
endif

if (Nbasics.gt.0) then
   do i=1,Nbasics
      Kb(i) = 10**(-pKb(i))
      fBHplus_bulk(i) =  Kb(i)/cOHmin / (1.0 + Kb(i)/cOHmin)
   enddo
endif

!! number density of charged beads in bulk !!

rhoacidsbulk = 0.

rhobasicsbulk = 0.

chargebalance = (xHplusbulk - xOHminbulk)/vsol 

do NC=1,Ncomp
  if (flagGC(NC).eq.1) then

    do i=1,long(NC)
      ia = acidtype(i,NC)
      ib = basictype(i,NC)
      rhoacidsbulk(ia,NC) = rhoacidsbulk(ia,NC) + rhopolbulk(NC) * fAmin_bulk(ia)  ! charged acid segments in bulk
      rhobasicsbulk(ib,NC) = rhobasicsbulk(ib,NC) + rhopolbulk(NC) * fBHplus_bulk(ib)  ! charged basic segments in bulk
    enddo

    do i=1,Nacids
       chargebalance = chargebalance - rhoacidsbulk(i,NC)   
    enddo
    do i=1,Nbasics
        chargebalance = chargebalance + rhobasicsbulk(i,NC)
    enddo

  endif
enddo

if(chargebalance.gt.0) then  ! excess positive charge in bulk
  xposbulk= rhosalt*vsol*vpos
  xnegbulk= rhosalt*vsol*vneg + chargebalance * vsol * vneg ! NaCl + HCl + CompCl 
else                         ! excess negative charge in bulk  
  xposbulk= rhosalt*vsol*vpos - chargebalance * vsol * vpos ! NaCl+ NaOH + NaComp
  xnegbulk= rhosalt*vsol*vneg
endif

!!! bulk volume fraction of the solvent !!!

xsolbulk=1-xposbulk-xnegbulk-xHplusbulk-xOHminbulk

do NC=1,Ncomp
   if (flagGC(NC).eq.1) xsolbulk = xsolbulk - xpolbulk(NC)
enddo


if (Nacids.gt.0) then
  do i=1,Nacids
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Files for multiple runs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do NC = 1, Ncomp
write(lnqfile,'(A4,I2.2,A4)')'lnq.',NC,'.dat'
open(unit=1533+NC,file=lnqfile)
write(rogfile,'(A4,I2.2,A4)')'rog.',NC,'.dat'
open(unit=2533+NC,file=rogfile)
enddo

end


subroutine calc_expmupol

use globals
use mcharge
use layer
use volume
use bulk
use longs
use MPI
use pis
use mkai

implicit none

integer NC, i, is, j, js, ia, ib !! dummy indices
real*8 totalvolpol


! bulk pro calculation


xpota_bulk(0) = 1.
xpotb_bulk(0) = 1.

do ia = 1,Nacids
  xpota_bulk(ia) = 1./fAmin_bulk(ia)
enddo

do ib = 1,Nbasics
  xpotb_bulk(ib) = 1./fBHplus_bulk(ib)
enddo

do is = 0, Npoorsv
  xpot_bulk(is) = xsolbulk**vpol(is)
  do js = 0, Npoorsv
    xpot_bulk(is) = xpot_bulk(is) * exp(rhopolbulk(NC) * sumaXu(is,js) * st(is,js))
  enddo
enddo

! expmupol calculation
expmupol = 0.0
probulk = 1.0

do NC=1,Ncomp
  if (flagGC(NC).eq.1) then
     expmupol(NC) = rhopolbulk(NC)*vsol
    
     do i=1,long(NC)

       is = segpoorsv(i,NC)
       ia = acidtype(i,NC)
       ib = basictype(i,NC)

       probulk(NC) = probulk(NC) / xpota_bulk(ia)
       probulk(NC) = probulk(NC) / xpotb_bulk(ib)   !! fraction of charged beads term of expmupol
       probulk(NC) = probulk(NC) / xpot_bulk(is)

       expmupol(NC) = expmupol(NC) / xpota_bulk(ia)
       expmupol(NC) = expmupol(NC) / xpotb_bulk(ib)   !! fraction of charged beads term of expmupol
       expmupol(NC) = expmupol(NC) / xpot_bulk(is)
        
     enddo
     expmupol(NC) = expmupol(NC) / totalcuantas(NC)    
     qbulk(NC) = probulk(NC) * totalcuantas(NC)
     sumprolnpro_bulk(NC) = probulk(NC)*dlog(probulk(NC)) 
     print*,expmupol(NC), qbulk(NC), sumprolnpro_bulk(NC)
  endif
   
enddo

end
