subroutine allocation
use globals
use mkinsol
use longs
use partfunc
use mkai

allocate (q(ntot))
allocate (avpol(0:Npoorsv, ntot))
allocate (xpol(ntot))
allocate (inn(0:Npoorsv,cuantas,ntot,base))
allocate (maxpos(cuantas,2*ntot))
allocate (minpos(cuantas,2*ntot))
allocate (eps(ntot))
allocate (xtotal(Npoorsv,2*ntot))
allocate (Xu(ntot,ntot))
allocate (pp(2*ntot))
end
