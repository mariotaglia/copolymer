subroutine read
use longs
use globals
use bulk
use MPI
implicit none
integer i

!     reading in of variables from stdin

read(8,*),nada
read(8,*),curvature

read(8,*),nada
read(8,*),ntot, maxntot

read(8,*),nada
read(8,*)cuantas

read(8,*),nada
read(8,*)long

READ(8,*),nada
read(8,*),nst
do i =1, nst
read(8,*),sts(i)
end do

READ(8,*),nada
read(8,*),infile

READ(8,*),nada
READ(8,*), nnpol

do i =1, nnpol
read(8,*),npols(i)
end do

read(8,*)nada
read(8,*)eps1

read(8,*)nada
read(8, *)Xulimit

read(8,*)nada
read(8,*)lseg

end
