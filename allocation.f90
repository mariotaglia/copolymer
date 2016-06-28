subroutine allocation
use globals
use mkinsol
use longs
use partfunc

allocate (q(ntot))
allocate (avpol(ntot,2))
allocate (xpol(ntot))
allocate (in1n(cuantas,ntot,base))
allocate (in2n(cuantas,ntot,base))
allocate (maxpos(cuantas,2*ntot))
allocate (minpos(cuantas,2*ntot))
allocate (eps(ntot))
allocate (xtotal(2*ntot))
allocate (Xu(ntot,ntot))
allocate (pp(2*ntot))
end
