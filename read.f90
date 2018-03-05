subroutine read
use mcharge
use longs
use globals
use bulk
use MPI
use mkai
implicit none
integer i, j
integer block_cuantas, restcuantas

!     reading in of variables from stdin

read(8,*)nada
read(8,*)curvature

read(8,*)nada
read(8,*)ntot, maxntotcounter_ini, maxntot

read(8,*)nada
read(8,*)totalcuantas

block_cuantas=int(totalcuantas/size/12)
cuantas=block_cuantas*12
restcuantas=totalcuantas-size*12*block_cuantas
if (rank.eq.(size-1))cuantas=cuantas+restcuantas

read(8,*)nada
read(8,*)long 

READ(8,*)nada
read(8,*)Npoorsv

read(8,*)nada
allocate(st(0:Npoorsv,0:Npoorsv))
st(0,0)=0.
do i = 1, Npoorsv
st(0,i)=0.
st(i,0)=0.
read(8,*)(st(i,j), j = 1, i)
 do j = 1, i
   st(j,i) = st(i,j)
 enddo
enddo

read(8,*)nada
allocate(dimf(0:Npoorsv,0:Npoorsv))
dimf(0,0)=0
do i = 1, Npoorsv
dimf(0,i)=0.
dimf(i,0)=0.
read(8,*)(dimf(i,j), j = 1, i)
 do j = 1, i
   dimf(j,i) = dimf(i,j)
 enddo
enddo

read(8,*)nada
read(8,*)Ncharge

read(8,*)nada
allocate(charge(Ncharge))
do i=1,Ncharge
read(8,*)charge(i)
enddo

READ(8,*)nada
read(8,*)infile

read(8,*)nada
read(8,*)flagkai

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

read(8,*)nada
read(8,*)Csalt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read chain structure from structure.in

open(file='structure.in', unit = 9)
allocate(segpoorsv(long))
allocate(chargetype(long))
do i = 1, long
read(9,*),segpoorsv(i), chargetype(i)
enddo

electroflag=0

do i=1,long
   if (chargetype(i).ne.0)electroflag=1
enddo
print*, "electroflag is", electroflag
end
