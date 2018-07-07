subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche



allocate (phi(0:2*ntot))
allocate (fAmin(Nacids,ntot), fBHplus(Nbasics,ntot), avpos(ntot), avneg(ntot), avHplus(ntot),avOHmin(ntot), xcharge(ntot))
allocate (Ntrans(long,cuantas))
allocate (trans(long))
allocate (Uchain(cuantas))
allocate (q(ntot))
allocate (sumprolnpro(ntot))
allocate (sumprouchain(ntot))
allocate (avpol(0:Npoorsv, ntot))
allocate (avpola(1:Nacids,ntot), avpolb(1:Nbasics,ntot))
allocate (xpol(ntot))
allocate (inn(0:Npoorsv,cuantas,ntot,base))
allocate (inn_a(0:Nacids,cuantas,ntot,base), inn_b(0:Nbasics,cuantas,ntot,base))
allocate (maxpos(cuantas,2*ntot))
allocate (minpos(cuantas,2*ntot))
allocate (xtotal(Npoorsv,2*ntot))
allocate (Xu(ntot,ntot,Npoorsv, Npoorsv))
allocate (xsol(ntot))
allocate (pp((npoorsv+2)*ntot))

ALLOCATE (epsfcn(0:ntot+1))
ALLOCATE (dielpol(1:ntot))
ALLOCATE (Depsfcn(0:ntot+1))

end
