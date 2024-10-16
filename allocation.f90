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
allocate (phi(0:dimR*2,dimZ))
allocate (fAmin(Nacids,dimR*2,dimZ), fBHplus(Nbasics,dimR*2,dimZ), avpos(dimR*2,dimZ), avneg(dimR*2,dimZ))
allocate (avHplus(dimR*2,dimZ),avOHmin(dimR*2,dimZ), xcharge(dimR*2,dimZ))
allocate (Ntrans(maxlong,cuantas_max,Ncomp))
allocate (trans(maxlong,Ncomp))
allocate (Uchain(cuantas_max, Ncomp))
allocate (q(dimR,dimZ,NComp))
allocate (sumprolnpro(dimR,dimZ,Ncomp))
allocate (sumprouchain(dimR,dimZ,Ncomp))
allocate (avpol(0:Npoorsv, dimR*2,dimZ,NComp))
allocate (avpola(0:Nacids,dimR*2,dimZ,Ncomp), avpolb(0:Nbasics,dimR*2,dimZ,Ncomp))
allocate (xpol(dimR*2,dimZ,Ncomp))
allocate (innZ(maxlong,cuantas_max,Ncomp),innR(maxlong,cuantas_max,2*maxntotR_max,Ncomp))
allocate (xtotal(Npoorsv,dimR*2,dimZ))

allocate (expmupol(Ncomp), rhopolbulk(Ncomp), xpolbulk(Ncomp), avpolbulk(0:Npoorsv,Ncomp))
allocate (probulk(Ncomp), sumprolnpro_bulk(Ncomp), qbulk(Ncomp))
allocate (fAmin_bulk(0:Nacids), fBHplus_bulk(0:Nbasics), rhoacidsbulk(0:Nacids,Ncomp), rhobasicsbulk(0:Nbasics,Ncomp))

allocate (xpot_bulk(0:Npoorsv),xpota_bulk(0:Nacids),xpotb_bulk(0:Nbasics))

Rini_kais=minntotR_min-10
Rfin_kais=maxntotR_max+10
if(Rini_kais.lt.1)Rini_kais = 1
if(Rfin_kais.gt.dimR)Rfin_kais = dimR
do NC=1,NComp
  if(flagGC(NC).eq.1) Rfin_kais = dimR + 100
enddo
allocate (Xu(1:Rfin_kais,1:Rfin_kais,-Xulimit:Xulimit,Npoorsv, Npoorsv), sumaXu(0:Npoorsv,0:Npoorsv), gtot(0:Npoorsv,0:Npoorsv))
allocate (xsol(dimR,dimZ))
allocate (pp((npoorsv+2)*ntot))
ALLOCATE (epsfcn(0:dimR+1,dimZ))
ALLOCATE (dielpol(dimR*2,dimZ))
ALLOCATE (Depsfcn(0:dimR*2,dimZ))
allocate (vchain(Ncomp))
end
