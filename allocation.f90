subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche
use volume
use bulk

allocate (Xflag((npoorsv+2)*ntot))
allocate (phi(0:dimR+dimbulk,dimZ))
allocate (fAmin(Nacids,dimR+dimbulk,dimZ), fBHplus(Nbasics,dimR+dimbulk,dimZ), avpos(dimR+dimbulk,dimZ), avneg(dimR+dimbulk,dimZ))
allocate (avHplus(dimR+dimbulk,dimZ),avOHmin(dimR+dimbulk,dimZ), xcharge(dimR+dimbulk,dimZ))
allocate (Ntrans(maxlong,cuantas_max,Ncomp))
allocate (trans(maxlong,Ncomp))
allocate (Uchain(cuantas_max, Ncomp))
allocate (q(dimR,dimZ,NComp))
allocate (sumprolnpro(dimR,dimZ,Ncomp))
allocate (sumprouchain(dimR,dimZ,Ncomp))
allocate (avpol(0:Npoorsv, dimR, dimZ, NComp))
allocate (avpola(0:Nacids, dimR, dimZ, Ncomp), avpolb(0:Nbasics, dimR, dimZ, Ncomp))
allocate (xpol(dimR,dimZ,Ncomp))
allocate (innZ(maxlong,cuantas_max,Ncomp),innR(maxlong,cuantas_max,2*maxntotR_max,Ncomp))
allocate (xtotal(Npoorsv,dimR+dimbulk,dimZ))

allocate (expmupol(Ncomp), rhopolbulk(Ncomp), xpolbulk(Ncomp), avpolbulk(0:Npoorsv,Ncomp))
allocate (probulk(Ncomp), sumprolnpro_bulk(Ncomp), qbulk(Ncomp))
allocate (fAmin_bulk(0:Nacids), fBHplus_bulk(0:Nbasics), rhoacidsbulk(0:Nacids,Ncomp), rhobasicsbulk(0:Nbasics,Ncomp))

allocate (xpot_bulk(0:Npoorsv),xpota_bulk(0:Nacids),xpotb_bulk(0:Nbasics))

Rini_kais=minntotR_min-10
Rfin_kais=maxntotR_max+10
if(Rini_kais.lt.1)Rini_kais = 1
if(Rfin_kais.gt.dimR)Rfin_kais = dimR
do NC=1,NComp
  if(flagGC(NC).eq.1) Rfin_kais = dimR + dimbulk
enddo
allocate (Xu(1:Rfin_kais,1:Rfin_kais,-Xulimit:Xulimit,Npoorsv, Npoorsv), sumaXu(0:Npoorsv,0:Npoorsv), gtot(0:Npoorsv,0:Npoorsv))
allocate (xsol(dimR,dimZ))
allocate (pp((npoorsv+2)*ntot))
ALLOCATE (epsfcn(0:dimR+1,dimZ))
ALLOCATE (dielpol(dimR+dimbulk,dimZ))
ALLOCATE (Depsfcn(0:dimR+dimbulk,dimZ))
allocate (vchain(Ncomp))
end
