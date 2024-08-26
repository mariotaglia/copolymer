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
allocate (phi(0:dimR+1,dimZ))
allocate (fAmin(Nacids,dimR,dimZ), fBHplus(Nbasics,dimR,dimZ), avpos(dimR,dimZ), avneg(dimR,dimZ))
allocate (avHplus(dimR,dimZ),avOHmin(dimR,dimZ), xcharge(dimR,dimZ))
allocate (Ntrans(maxlong,cuantas_max,Ncomp))
allocate (trans(maxlong,Ncomp))
allocate (Uchain(cuantas_max, Ncomp))
allocate (q(dimR,dimZ,NComp))
allocate (sumprolnpro(dimR,dimZ,Ncomp))
allocate (sumprouchain(dimR,dimZ,Ncomp))
allocate (avpol(0:Npoorsv, dimR,dimZ,NComp))
allocate (avpola(0:Nacids,dimR,dimZ,Ncomp), avpolb(0:Nbasics,dimR,dimZ,Ncomp))
allocate (xpol(dimR,dimZ,Ncomp))
allocate (innZ(maxlong,cuantas_max,Ncomp),innR(maxlong,cuantas_max,maxntotR_max,Ncomp))
allocate (xtotal(Npoorsv,dimR,dimZ))

allocate (expmupol(Ncomp), rhopolbulk(Ncomp), xpolbulk(Ncomp), avpolbulk(0:Npoorsv,Ncomp))
allocate (fAmin_bulk(0:Nacids), fBHplus_bulk(0:Nbasics), rhoacidsbulk(0:Nacids,Ncomp), rhobasicsbulk(0:Nbasics,Ncomp))

allocate (xpot_bulk(0:Npoorsv),xpota_bulk(0:Nacids),xpotb_bulk(0:Nbasics))

Rini_kais=minntotR_min-10
Rfin_kais=maxntotR_max+10
if(Rini_kais.lt.1)Rini_kais = 1
if(Rfin_kais.gt.dimR)Rfin_kais = dimR
allocate (Xu(dimR,Rini_kais:Rfin_kais,-Xulimit:Xulimit,Npoorsv, Npoorsv), sumaXu(0:Npoorsv,0:Npoorsv))
allocate (xsol(dimR,dimZ))
allocate (pp((npoorsv+2)*ntot))
ALLOCATE (epsfcn(0:dimR+1,dimZ))
ALLOCATE (dielpol(dimR,dimZ))
ALLOCATE (Depsfcn(0:dimR+1,dimZ))
allocate (vchain(Ncomp))
end
