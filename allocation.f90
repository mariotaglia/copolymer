subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche



allocate (phi(0:dimR,0:dimZ))
allocate (fAmin(Nacids,dimR,dimZ), fBHplus(Nbasics,dimR,dimZ), avpos(dimR,dimZ), avneg(dimR,dimZ))
allocate (avHplus(dimR,dimZ),avOHmin(dimR,dimZ), xcharge(dimR,dimZ))
allocate (Ntrans(long,cuantas))
allocate (trans(long))
allocate (Uchain(cuantas))
allocate (q(dimR,dimZ))
allocate (sumprolnpro(dimR,dimZ))
allocate (sumprouchain(dimR,dimZ))
allocate (avpol(0:Npoorsv, dimR,dimZ))
allocate (avpola(1:Nacids,dimR,dimZ), avpolb(1:Nbasics,dimR,dimZ))
allocate (xpol(ntot))
allocate (inn(0:Npoorsv,cuantas,ntot,base))
allocate (inn_a(0:Nacids,cuantas,ntot,base), inn_b(0:Nbasics,cuantas,ntot,base))
allocate (maxpos(cuantas,2*ntot))
allocate (minpos(cuantas,2*ntot))
allocate (xtotal(Npoorsv,ntot))
allocate (Xu(40,40,-Xulimit:Xulimit,Npoorsv, Npoorsv))
allocate (xsol(dimR,dimR))
allocate (pp((npoorsv+2)*ntot))

ALLOCATE (epsfcn(0:ntot+1))
ALLOCATE (dielpol(1:ntot))
ALLOCATE (Depsfcn(0:ntot+1))

end
