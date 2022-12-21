subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche
use volume

allocate (Xflag((npoorsv+2)*ntot))
allocate (phi(0:dimR+1,dimZ))
allocate (fAmin(Nacids,dimR,dimZ), fBHplus(Nbasics,dimR,dimZ), avpos(dimR,dimZ), avneg(dimR,dimZ))
allocate (avHplus(dimR,dimZ),avOHmin(dimR,dimZ), xcharge(dimR,dimZ))
allocate (Ntrans(maxlong,cuantas,Ncomp))
allocate (trans(maxlong,Ncomp))
allocate (Uchain(cuantas, Ncomp))
allocate (q(dimR,dimZ,NComp))
allocate (sumprolnpro(dimR,dimZ,Ncomp))
allocate (sumprouchain(dimR,dimZ,Ncomp))
allocate (avpol(0:Npoorsv, dimR,dimZ,NComp))
allocate (avpola(0:Nacids,dimR,dimZ,Ncomp), avpolb(0:Nbasics,dimR,dimZ,Ncomp))
allocate (xpol(dimR,dimZ,Ncomp))
allocate (innZ(maxlong,cuantas,Ncomp),innR(maxlong,cuantas,maxntotR,Ncomp))
allocate (xtotal(Npoorsv,dimR,dimZ))

Rini_kais=minntotR-10
Rfin_kais=maxntotR+10
if(Rini_kais.lt.1)Rini_kais = 1
if(Rfin_kais.gt.dimR)Rini_kais = dimR

allocate (Xu(Rini_kais:Rfin_kais,Rini_kais:Rfin_kais,-Xulimit:Xulimit,Npoorsv, Npoorsv))
allocate (xsol(dimR,dimZ))
allocate (pp((npoorsv+2)*ntot))
ALLOCATE (epsfcn(0:dimR+1,dimZ))
ALLOCATE (dielpol(dimR,dimZ))
ALLOCATE (Depsfcn(0:dimR+1,dimZ))
allocate (vchain(Ncomp))
end
