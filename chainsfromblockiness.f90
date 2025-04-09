subroutine genchains(seed,lambda, NA, NB) 

use blockiness

implicit none

real*8 lambda
real*8 rands, generator
integer seed
integer NA, NB, Ntotal
real*8 fA, fB, pAA, pAB, pBB, pBA

generator = rands(seed)
Ntotal = NA + NB

fA = NA/Ntotal
print*,fA

pAA = fA*(1. - lambda) + lambda
pBB = fA*(lambda - 1.) + 1.
pBA = 1. - pAA
pAB = 1. - pBB


