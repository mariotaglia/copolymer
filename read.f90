subroutine read
use mcharge
use longs
use globals
use bulk
use MPI
use mkai
use volume

implicit none
integer i, j
integer block_cuantas, restcuantas

!     reading in of variables from stdin


read(8,*)nada
read(8,*)curvature

read(8,*)nada
read(8,*)dimR, dimZ, maxntotR, maxntotZ

ntot=dimR*dimZ

read(8,*)nada
read(8,*)totalcuantas

block_cuantas=int(totalcuantas/size/12)
cuantas=block_cuantas*12
restcuantas=totalcuantas-size*12*block_cuantas
if (rank.eq.(size-1))cuantas=cuantas+restcuantas

read(8,*)nada
read(8,*)long !number of segments 

READ(8,*)nada
read(8,*)Npoorsv !number of types of poor sv segments

allocate(vpol(0:Npoorsv))

read(8,*)nada
do i=0, Npoorsv
  read(8,*)vpol(i)
enddo
!read(8,*)nada
!allocate(st(0:Npoorsv,0:Npoorsv))
!st(0,0)=0.
!do i = 1, Npoorsv
!st(0,i)=0.
!st(i,0)=0.
!read(8,*)(st(i,j), j = 1, i)
! do j = 1, i
!   st(j,i) = st(i,j)
! enddo
!enddo

read(8,*)nada
allocate(dimfkais(0:Npoorsv,0:Npoorsv),dimf(0:Npoorsv,0:Npoorsv))
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
read(8,*)Nacids, Nbasics ! number of types of acid segments and basic segments

allocate(vpol_a(Nacids), vpol_b(Nbasics),Ka(Nacids), Kb(Nbasics), pKa(Nacids), pKb(Nbasics))

pKa=0.0
pKb=0.0
vpol_a=0.0
vpol_b=0.0

read(8,*)nada
do i=1,Nacids
  read(8,*)vpol_a(i)
enddo

read(8,*)nada
do i=1,Nbasics
  read(8,*)vpol_b(i)
enddo

if (Nacids.gt.0) then
  read(8,*)nada
  do i=1,Nacids
    read(8,*)pKa(i) ! acid constants of each acid segment
  enddo
endif

if (Nbasics.gt.0) then
  read(8,*)nada
  do i=1,Nbasics ! basic constants of each basic segment
    read(8,*)pKb(i)
  enddo
endif

READ(8,*)nada
read(8,*)infile

read(8,*)nada
read(8,*)flagkai

READ(8,*)nada
READ(8,*)npolini, npolfirst, npollast, npolstep

read(8,*)nada
read(8,*)Xulimit

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

read(8,*)nada
read(8,*)pHbulk

read(8,*)nada
read(8,*)dielP

read(8,*)nada
read(8,*)nbranches ! number of branches

allocate (branch_pos(nbranches))
allocate (branch_long(nbranches))
long_branches = 0

do j = 1, nbranches
read(8,*) branch_pos(j), branch_long(j)
long_branches = long_branches + branch_long(j)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read chain structure from structure.in

open(file='structure.in', unit = 9)
allocate(segpoorsv(long))
allocate(acidtype(long))
allocate(basictype(long))
allocate(torsionstate(long))
do i = 1, long
read(9,*),segpoorsv(i), acidtype(i), basictype(i), torsionstate(i) 
print torsionstate(i)
enddo
stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read st(i,j) from epsilon.in

open(file='epsilon.in', unit=10)
allocate(st(0:Npoorsv,0:Npoorsv))
st(0,0)=0.
do i = 1, Npoorsv
st(0,i)=0.
st(i,0)=0.
read(10,*)(st(i,j), j = 1, i)
 do j = 1, i
   st(j,i) = st(i,j)
 enddo
enddo

end
