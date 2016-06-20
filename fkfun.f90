subroutine fkfun(x,f,ier2)
use globals
use partfunc
use layer
use volume
use bulk
use longs
use MPI
use pis
implicit none

integer*4 ier2
real*8 protemp
real*8 x(2*ntot),f(2*ntot)
real*8 xh(2*ntot)
real*8 xpot(2*ntot,2)
real*8 pro(cuantas)
integer k,i,j,k1,k2,ii, jj,iz       ! dummy indices
integer err
integer n
real*8 avpol_tmp(2*ntot,2)
real*8 avpol_tosend(ntot,2)
real*8 algo, algo1,algo2
double precision, external :: factorcurv
real*8 sumpol

! Jefe
if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, 2*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
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


do i = 1, ntot
protemp = dlog(xh(i)**(vpol))
xpot(i,1) = dexp(protemp)
  do j = 1, ntot
      protemp = protemp+st*vsol*vpol*Xu(i,j)*xtotal(j)
  end do
      protemp = protemp+eps(i)
xpot(i,2) = dexp(protemp)
enddo

xpot(n+1:2*n,1)=xpot(n,1)
xpot(n+1:2*n,2)=xpot(n,2)

!    probability distribution
q=0.0d0                   ! init q to zero

do ii=1,ntot ! position of segment #0 
 do i=1,cuantas

 pro(i) = 1.0

    do j=minpos(i,ii), maxpos(i,ii) ! posicion dentro del poro
     k = j-minpos(i,ii)+1 ! k may be lager than ntot
     pro(i)= pro(i) * xpot(j,1)**in1n(i,ii,k) ! hidrofilico
     pro(i)= pro(i) * xpot(j,2)**in2n(i,ii,k)
    enddo

    q=q+pro(i)

    do j=minpos(i,ii), maxpos(i,ii)
     k = j-minpos(i,ii)+1 ! k may be larger than ntot
     avpol_tmp(j,1)=avpol_tmp(j,1)+pro(i)*vpol*in1n(i,ii,k)*factorcurv(ii,j)
     avpol_tmp(j,2)=avpol_tmp(j,2)+pro(i)*vpol*in2n(i,ii,k)*factorcurv(ii,j)
    enddo

 enddo ! i
enddo   ! ii

avpol_tosend(1:ntot,:)=avpol_tmp(1:ntot,:) 

!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, 2*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

! norma avpol
! integrate over whole system

sumpol = 0.0

do i = 1, ntot
select case (curvature)
case (0)
sumpol = sumpol + (avpol(i,1) + avpol(i,2))*delta ! final result in units of chains/nm^2
case(1)
sumpol = sumpol + (avpol(i,1) + avpol(i,2))*(float(i)-0.5)*delta*delta*2.0*pi ! final result in units of chains/nm
case(2)
sumpol = sumpol + (avpol(i,1) + avpol(i,2))*(((float(i)-0.5)*delta)**2)*delta*4.0*pi ! final result in units of chains/micelle
end select
enddo

sumpol = sumpol/(vpol*vsol)/long
avpol = avpol/sumpol*npol ! integral of avpol is fixed

! contruction of f and the volume fractions

do i=1,n
 f(i)=xh(i)-1.0d0
 f(i) = f(i) + avpol(i,1)+avpol(i,2)
enddo

do i = 1,n ! xtotal
! f(i+n) = 0.0
 f(i+n) = -avpol(i,2)+xtotal(i)
enddo

iter=iter+1

algo = 0.0
algo1 = 0.0
algo2 = 0.0
do i = 1, n*2
! algo1 = algo1 + f(i)**2
! algo2 = algo2 + f(i+n)**2
algo = algo + f(i)**2
end do

if(rank.eq.0)PRINT*, iter, algo, sumpol
!if(rank.eq.0)PRINT*, iter, algo1,algo2,algo1+algo2
norma=algo

3333 continue

return
end
