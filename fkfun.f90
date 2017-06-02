subroutine fkfun(x,f,ier2)
use globals
use partfunc
use layer
use volume
use bulk
use longs
use MPI
use pis
use mkai
implicit none

integer*4 ier2
real*8 protemp, sttemp
real*8 x((Npoorsv+1)*ntot),f((Npoorsv+1)*ntot)
real*8 xh(2*ntot)
real*8 xpot(0:Npoorsv,2*ntot)
real*8 pro(cuantas)
integer k,i,j,k1,k2,ii, jj,iz       ! dummy indices
integer is, js
integer err
integer n
real*8 avpol_tmp(0:Npoorsv,2*ntot)
real*8 avpol_tosend(0:Npoorsv, ntot)
real*8 xpol_tosend(ntot)
real*8 algo, algo1,algo2
double precision, external :: factorcurv
real*8 sumpol
real*8 q_tosend(ntot)
real*8 qall_tosend

! Jefe
if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, (Npoorsv+1)*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

n = ntot

! chain parameters

do i=1,n                 
xh(i)=x(i)
enddo

do i=1,n
do is = 1,Npoorsv                
xtotal(is,i) = x(i+n*is) 
enddo
enddo

do i = n+1,2*n
xtotal(:,i) = 0.0 ! bulk
xh(i) = 1.0
enddo

!sttemp = st/(vpol*vsol)

do i = 1, ntot
protemp = dlog(xh(i)**(vpol))
xpot(0,i) = dexp(protemp)
enddo 



do is = 1, Npoorsv
  do i = 1, ntot

! calculate xpot(i, is)

  protemp = 0.0

   do js = 1, Npoorsv 
   do j = 1, ntot
      protemp = protemp+st(is,js)/(vpol*vsol)*Xu(i,j)*xtotal(js,j)
   enddo
   enddo

   xpot(is,i) = xpot(0,i)*dexp(protemp)
 enddo
enddo

do is = 0,Npoorsv
xpot(is,n+1:2*n)=xpot(is,n)
enddo

!    probability distribution

avpol_tosend = 0.0
xpol_tosend = 0.0
avpol_tmp = 0.0
avpol = 0.0
xpol = 0.0
q = 0.0
q_tosend=0.0d0                   ! init q to zero
qall_tosend = 0.0
qall = 0.0


  do ii=1,maxntotcounter ! position of center of mass 
   do i=1,cuantas ! loop over conformations

     pro(i) = 1.0

     do j=minpos(i,ii), maxpos(i,ii) ! loop over lattice position
      k = j-minpos(i,ii)+1 ! k may be lager than ntot

      do is = 0, Npoorsv 
      pro(i)= pro(i) * xpot(is,j)**inn(is,i,ii,k) 
      enddo
     enddo

      q_tosend(ii)=q_tosend(ii)+pro(i)
      qall_tosend = qall_tosend + pro(i)

      xpol_tosend(ii)=xpol_tosend(ii)+pro(i)

     do j=minpos(i,ii), maxpos(i,ii)
      k = j-minpos(i,ii)+1 ! k may be larger than ntot

      do is = 0, Npoorsv 
      avpol_tmp(is,j)=avpol_tmp(is,j)+pro(i)*vpol*inn(is,i,ii,k)*factorcurv(ii,j)
      enddo

     enddo

   enddo ! i
  enddo   ! ii


avpol_tosend(:, 1:ntot)=avpol_tmp(:, 1:ntot) 

!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, (Npoorsv+1)*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(xpol_tosend, xpol, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(q_tosend, q, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(qall_tosend, qall, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, (Npoorsv+1)*ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(xpol_tosend, xpol, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(q_tosend, q, ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  call MPI_REDUCE(qall_tosend, qall, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  ier2 = 0
  goto 3333
endif


! norma avpol
! integrate over whole system

sumpol = 0.0

do i = 1, ntot
do is = 0, Npoorsv

select case (curvature)
case (0)
sumpol = sumpol + avpol(is,i)*delta ! final result in units of chains/nm^2
case(1)
sumpol = sumpol + avpol(is,i)*(float(i)-0.5)*delta*delta*2.0*pi ! final result in units of chains/nm
case(2)
sumpol = sumpol + avpol(is,i)*(((float(i)-0.5)*delta)**2)*delta*4.0*pi ! final result in units of chains/micelle
end select

enddo
enddo


sumpol = sumpol/(vpol*vsol)/long
avpol = avpol/sumpol*npol ! integral of avpol is fixed

sumpol = 0.0
do i = 1, ntot
select case (curvature)
case (0)
sumpol = sumpol + xpol(i)*delta ! final result in units of chains/nm^2
case(1)
sumpol = sumpol + xpol(i)*(float(i)-0.5)*delta*delta*2.0*pi ! final result in units of chains/nm
case(2)
sumpol = sumpol + xpol(i)*(((float(i)-0.5)*delta)**2)*delta*4.0*pi ! final result in units of chains/micelle
end select
enddo
xpol = xpol/sumpol*npol ! integral of avpol is fixed




! contruction of f and the volume fractions

do i=1,n
 f(i)=xh(i)-1.0d0
 do is=0, Npoorsv
   f(i) = f(i) + avpol(is,i)
 enddo
enddo

do is=1,Npoorsv
do i=1,n ! xtotal
 f(i+n*is) = -avpol(is,i)+xtotal(is,i)
enddo
enddo

iter=iter+1

algo = 0.0
algo1 = 0.0
algo2 = 0.0
do i = 1, n*(Npoorsv+1)
algo = algo + f(i)**2
end do

if(rank.eq.0)PRINT*, iter, algo, q(1), Q(2), q(3), q(4)
norma=algo

3333 continue


ier2 = 0

return
end
