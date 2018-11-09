subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche
use volume


allocate (phi(0:dimR+1,dimZ))
allocate (fAmin(Nacids,dimR,dimZ), fBHplus(Nbasics,dimR,dimZ), avpos(dimR,dimZ), avneg(dimR,dimZ))
allocate (avHplus(dimR,dimZ),avOHmin(dimR,dimZ), xcharge(dimR,dimZ))
allocate (Ntrans(long,cuantas,Ncomp))
allocate (trans(long))
allocate (Uchain(cuantas, Ncomp))
allocate (q(dimR,dimZ,NComp))
allocate (sumprolnpro(dimR,dimZ,Ncomp))
allocate (sumprouchain(dimR,dimZ,Ncomp))
allocate (avpol(0:Npoorsv, dimR,dimZ,NComp))
allocate (avpola(0:Nacids,dimR,dimZ,Ncomp), avpolb(0:Nbasics,dimR,dimZ,Ncomp))
allocate (xpol(dimR,dimZ))
allocate (innZ(long,cuantas),innR(long,cuantas,maxntotR))
allocate (xtotal(Npoorsv,dimR,dimZ))
allocate (Xu(dimR,dimR,-Xulimit:Xulimit,Npoorsv, Npoorsv))
allocate (xsol(dimR,dimR))
allocate (pp((npoorsv+2)*ntot))
ALLOCATE (epsfcn(0:dimR+1,dimZ))
ALLOCATE (dielpol(dimR,dimZ))
ALLOCATE (Depsfcn(0:dimR+1,dimZ))
allocate (vchain(Ncomp))
end
