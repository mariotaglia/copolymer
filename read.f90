subroutine read
use longs
use globals
use bulk
use MPI
use mkai
implicit none
integer i, j

!     reading in of variables from stdin

read(8,*)nada
read(8,*)curvature

read(8,*)nada
read(8,*)ntot, maxntotcounter_ini, maxntot

read(8,*)nada
read(8,*)cuantas

read(8,*)nada
read(8,*)long 

READ(8,*)nada
read(8,*)Npoorsv

read(8,*)nada
allocate(st(Npoorsv,Npoorsv))
do i = 1, Npoorsv
read(8,*)(st(i,j), j = 1, i)
 do j = 1, i
   st(j,i) = st(i,j)
 enddo
enddo

READ(8,*)nada
read(8,*)infile

READ(8,*)nada
READ(8,*)npolini, npolfirst, npollast, npolstep

read(8,*)nada
read(8,*)eps1

read(8,*)nada
read(8, *)Xulimit

read(8,*)nada
read(8,*)lseg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read chain structure from structure.dat

open(file='structure.dat', unit = 9)
allocate(segpoorsv(long))

do i = 1, long
read(9,*),segpoorsv(i)
enddo
end
