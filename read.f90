subroutine read
use longs
use globals
use bulk
use MPI
use kai
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
read(8,*)long, long1

READ(8,*)nada
read(8,*)Npoorsv

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

end
