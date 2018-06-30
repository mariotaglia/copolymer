subroutine print_ent2(xend)

use longs
implicit none

real*8 xend(3,200)
integer i,j, k
character*21 filename
integer, save :: indexncha = 1

! Imprime cadenas en formato ENT

write(filename,'(A6,A1, I3.3, A4)') 'cadena','.', indexncha, '.ent'

open(unit=4400, file=filename)

do i=1, long ! Imprime todo
WRITE(4400,'(A6,I5,A3,I12,A4,F8.3,F8.3,F8.3)') &
"HETATM",i,"  C",i,"    ",xend(1, i)*10,  &
xend(2, i)*10,xend(3, i)*10
end do

   WRITE((4400),'(A6,I5,I5)')"CONECT", 1, 2

i = 2
do j = 3, long-long_branches ! Une segmentos
  i = i + 1
  WRITE((4400),'(A6,I5,I5)')"CONECT", i-1, i
end do

do k = 1, nbranches
 i = i + 1
 WRITE((4400),'(A6,I5,I5)')"CONECT", branch_pos(k), i

 do j = 2, branch_long(k)
   i = i + 1
   WRITE((4400),'(A6,I5,I5)')"CONECT", i-1, i
 enddo ! j
enddo ! k

WRITE(4400,*)"END"
close(4400)
indexncha = indexncha + 1
if(indexncha.eq.1000)stop
end subroutine

