subroutine allocation
use globals
use mcharge
use mkinsol
use longs
use partfunc
use mkai
use transgauche



allocate (phi(0:2*ntot))
allocate (avpos(ntot),avneg(ntot),xcharge(ntot))
!allocate (Ntrans(long,cuantas))
allocate (trans(long))
!allocate (Uchain(cuantas))
allocate (q(ntot))
allocate (sumprolnpro(ntot))
allocate (sumprouchain(ntot))
allocate (avpol(0:Npoorsv+2, ntot))
allocate (avpolc(1:Ncharge,ntot))
allocate (xpol(ntot))
!allocate (inn(0:Npoorsv,cuantas,ntot,base))
!allocate (innc(0:Ncharge,cuantas,ntot,base))
!allocate (maxpos(cuantas,2*ntot))
!allocate (minpos(cuantas,2*ntot))
allocate (eps(ntot))
allocate (xtotal(Npoorsv,2*ntot))
allocate (Xu(ntot,ntot,Npoorsv, Npoorsv))
allocate (xsol(ntot))
allocate (pp((npoorsv+2)*ntot))

end
