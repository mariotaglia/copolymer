subroutine fkfun(x,f,ier2)
use multicapa
use partfunc
use layer
use volume
use bulk
use longs
use MPI
implicit none

integer*4 ier2
real*8 protemp
real*8 x(ntot),f(ntot)
real*8 xh(2*ntot)
real*8 xpot(2*ntot,2)
real*8 pro(cuantas)
integer i,j,k1,k2,ii, jj,iz       ! dummy indices
integer err
integer n
real*8 avpol_tmp(ntot)

double precision, external :: factorcurv

! Jefe
if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

n = ntot

! chain parameters

do i=1,n                 
xh(i)=x(i)
xtotal(i) = x(i+n) 
enddo

do i = n+1,2*n
xtotal(i) = 0.0 ! bulk
xh(i) = 1.0
enddo

xpot(i,1) = dexp(protemp)

do i = 1, 2*ntot
protemp = dlog(xh(i)**(vpol))
  do j = 1, ntot
      protemp = protemp+st*vsol*vpol*Xu(i,j)*xtotal(j)
  end do
enddo
xpot(i,2) = dexp(protemp)
enddo


!    probability distribution
q=0.0d0                   ! init q to zero

do ii=1,ntot ! position of segment #0 
 do i=1,cuantas

 pro(i) = 1.0

    do j=1, ntot ! posicion
     pro(i)= pro(i) * xpot(j,1)**in1n(i,ii,j) ! hidrofilico
     pro(i)= pro(i) * xpot(j,2)**in2n(i,ii,j)
    enddo

    q=q+pro(i)

    do j=1,ntot
     avpol_tmp(j,1)=avpol_tmp(j,1)+pro(i)*vpol*in1n(i,ii,j)*factorcurv(ii,j)
     avpol_tmp(j,2)=avpol_tmp(j,2)+pro(i)*vpol*in2n(i,ii,j)*factorcurv(ii,j)
    enddo

 enddo ! i
enddo   ! ii

!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tmp, avpol, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tmp, avpol, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

! norma avpol

sumpol = 0.0
do i = 1, ntot
sumpol = sumpol + avpol(i,1) + avpol(i,2)
enddo
sumpol = sumpol/(vpol*vsol)/long
avpol = avpol/sumpol*nnpol ! integral of avpol is fixed


! contruction of f and the volume fractions

do i=1,n
 f(i)=xh(i)-1.0d0
 f(i) = f(i) + avpol(i,1)+avpol(i,2)
enddo

do i = 1,n ! xtotal
 f(i+n) = f(i+n) + avpol(i,2)-xtotal(i)
enddo

iter=iter+1

algo = 0.0
do i = 1, 2*n
 algo = algo + f(i)**2
end do

if(rank.eq.0)PRINT*, iter, algo
norma=algo

3333 continue

return
end
