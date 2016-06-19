subroutine allocation
use globals
use mkinsol
use longs

allocate (avpol(ntot,2))
allocate (in1n(maxcuantas,ntot,base))
allocate (in2n(maxcuantas,ntot,base))
allocate (maxpos(maxcuantas,2*ntot))
allocate (minpos(maxcuantas,2*ntot))
allocate (eps(ntot))
allocate (xtotal(2*ntot))
allocate (Xu(ntot,ntot))
allocate (pp(2*ntot))
end
