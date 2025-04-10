subroutine genchains(seed, lambda, beadcountA, beadcountB, beadtosend) 

use modblockiness

implicit none

real*8 lambda, lambdac, lambdae
real*8 rands, generator
integer seed, firstidx
integer beadcountA, beadcountB
integer nA, nB, nAT, nBT, nAA, nBA, nAB, nBB
real*8 fA, fB, pAA, pAB, pBB, pBA, wA, wB
real*8 pAAc, pABc, pBBc, pBAc
integer beadtosend(Ntotal)
integer i, idx ! dummy indexes

nAT = beadcountA
nBT = beadcountB

fA = float(beadcountA)/float(Ntotal)


pAA = fA*(1. - lambda) + lambda
pBB = fA*(lambda - 1.) + 1.
pBA = 1. - pAA
pAB = 1. - pBB

nA = 0
nB = 0
nAA = 0
nAB = 0
nBA = 0
nBB = 0

firstidx = 1

if (rands(seed).lt.fA) then
   nA = nA + 1
   nAT = nAT - 1
   beadtosend(firstidx) = 0

else
   nB = nB + 1
   nBT = nBT - 1
   beadtosend(firstidx) = 1
endif

do i = 2, Ntotal
  
  if (beadtosend(i - 1).eq.0) then
     wA = float(nAT) / float(nAT + nBT) 
     wA = wA / (pAA * float(nAT) / float(nAT + nBT) + (pAB * float(nBT) / float(nAT + nBT)))
     if (rands(seed).lt.(pAA * wA)) then
        beadtosend(i) = 0
        nA = nA + 1
        nAT = nAT - 1
        nAA = nAA + 1
     else
        beadtosend(i) = 1
        nB = nB + 1
        nBT = nBT - 1
        nAB = nAB + 1
     endif   
  elseif (beadtosend(i - 1).eq.1) then
     wB = float(nBT) / float(nAT + nBT) 
     wB = wB / (pBB * float(nBT) / float(nAT + nBT) + pBA * (float(nAT) / float(nAT + nBT)))
     if (rands(seed).lt.(Pbb * wB)) then
        beadtosend(i) = 1
        nB = nB + 1
        nBT = nBT - 1
        nBB = nBB + 1
     else
        beadtosend(i) = 0
        nA = nA + 1
        nAT = nAT -1
        nBA = nBA + 1
     endif
  endif
enddo

nAT = nAA + nAB
nBT = nBB + nBA

pAAc = float(nAA) / float(nA)
pBBc = float(nBB) / float(nB)
pABc = float(nAB) / float(nA)
pBAc = float(nBA) / float(nB)
lambdac = pAAc * pBBc - pABc * pBAc
lambdae = pAA * pBB - pBA * pAB

return
end
