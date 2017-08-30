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

read(8,*)nada
allocate(dimf(Npoorsv,Npoorsv))
do i = 1, Npoorsv
read(8,*)(dimf(i,j), j = 1, i)
 do j = 1, i
   dimf(j,i) = dimf(i,j)
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

read(8,*)nada
allocate(Ut(0:Npoorsv))
allocate(Ug(0:Npoorsv))
do i = 0, Npoorsv
read(8,*), Ut(i), Ug(i)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read chain structure from structure.in

open(file='structure.in', unit = 9)
allocate(segpoorsv(long))

do i = 1, long
read(9,*),segpoorsv(i)
enddo
end
