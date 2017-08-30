subroutine allocation
use globals
use mkinsol
use longs
use partfunc
use mkai
use transgauche



allocate (Ntrans(long,cuantas))
allocate (trans(long))
allocate (Uchain(cuantas))
allocate (q(ntot))
allocate (sumprolnpro(ntot))
allocate (sumprouchain(ntot))
allocate (avpol(0:Npoorsv, ntot))
allocate (xpol(ntot))
allocate (inn(0:Npoorsv,cuantas,ntot,base))
allocate (maxpos(cuantas,2*ntot))
allocate (minpos(cuantas,2*ntot))
allocate (eps(ntot))
allocate (xtotal(Npoorsv,2*ntot))
allocate (Xu(ntot,ntot,Npoorsv, Npoorsv))
allocate (xsol(ntot))
allocate (pp((npoorsv+1)*ntot))

end
