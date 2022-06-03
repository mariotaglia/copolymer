subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche
use volume

allocate (Xflag((npoorsv+4)*ntot))  !LEO
allocate (phi(0:dimR+1,dimZ))
!allocate (fAmin(Nacids,dimR,dimZ), fBHplus(Nbasics,dimR,dimZ), avpos(dimR,dimZ), avneg(dimR,dimZ)) !LEO
allocate (avpos(dimR,dimZ), avneg(dimR,dimZ)) ! LEO
allocate (avHplus(dimR,dimZ),avOHmin(dimR,dimZ), xcharge(dimR,dimZ))
allocate (fcopANC(dimR,dimZ), fcopAC(dimR,dimZ), fcopAion(dimR,dimZ), fASmol(dimR,dimZ),fAScopA(dimR,dimZ))!LEO
allocate (fmolNC(dimR,dimZ), fmolC(dimR,dimZ), fmolion(dimR,dimZ)) !LEO
allocate (xNcopA(dimR,dimZ), xNmol(dimR,dimZ)) ! LEO
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
allocate (Xu(dimR,dimR,-Xulimit:Xulimit,Npoorsv, Npoorsv))
allocate (xsol(dimR,dimZ))
allocate (pp((npoorsv+4)*ntot)) !LEO
ALLOCATE (epsfcn(0:dimR+1,dimZ))
ALLOCATE (dielpol(dimR,dimZ))
ALLOCATE (Depsfcn(0:dimR+1,dimZ))
allocate (vchain(Ncomp))
end
