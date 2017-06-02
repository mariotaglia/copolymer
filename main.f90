
call initmpi
call read
call allocation
call kai
call solve
end

subroutine solve
use globals
use partfunc
use layer
use volume
use bulk
use seed1
use longs
use MPI
use mkai

implicit none

integer is
integer *4 ier ! Kinsol error flag
real*8 pi
real*8 Na               
parameter (Na=6.02d23)
integer av1(ntot), av2(ntot)
real*8 avtmp
real*8 x1((npoorsv+1)*ntot),xg1((npoorsv+1)*ntot),x1ini((npoorsv+1)*ntot)   ! density solvent iteration vector
real*8 zc(ntot)           ! z-coordinate layer 

REAL*8 sumrhoz, meanz     ! Espesor medio pesado
real*8 pro                ! probability distribution function 
real*8 trash

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 fnorm              ! L2 norm of residual vector function fcn

external fcnelect         ! function containing the SCMFT eqs for solver
integer i,j,k,m,ii,flag,c, jj ! dummy indice0s

INTEGER temp
real*8 tempr
real*8 tmp

real*8 min1               ! variable to determine minimal position of chain
integer qqq,www,eee

integer il,inda,ncha

REAL*8 xfile((npoorsv+1)*ntot)                        
real*8 algo, algo2                  


integer*1 in1(long)
real*8 chains(3,long,ncha_max) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
real*8 zp(long)

real*8 sum,sumel          ! auxiliary variable used in free energy computation  
real*8 sumpi,sumrho,sumrhopol, sumrho2, sumrho2mol !suma de la fraccion de polimero

! global running files
!character*15 meanzfilename
!character*15 sigmafilename
!character*17 sigmaadfilename

! single layer files
character*18 sysfilename      ! contains value of free energy, input parameter etc
character*26 denssolfilename  ! contains the denisty of the solvent
character*27 lnqfilename  ! contains the denisty of the solvent
character*28 densendfilename
CHARACTER*24 totalfilename
CHARACTER*24 xtotalfilename
character*50, allocatable :: denspolfilename(:)

integer countfile         ! enumerates the outputfiles 
integer conf              ! counts number of conformations

integer readsalt          !integer to read salt concentrations


INTEGER cc

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

integer in1tmp(long)

allocate(denspolfilename(0:Npoorsv))

error = 1.0d-6
seed=435+ 3232*rank               ! seed for random number generator

print*, 'I am', rank, ' and my seed is', seed

if(rank.eq.0)print*, 'Program Multicapa'
if(rank.eq.0)print*, 'GIT Version: ', _VERSION

!     initializations of variables 
pi=dacos(-1.0d0)          ! pi = arccos(-1) 
itmax=200                 ! maximum number of iterations       
n=ntot                    ! size of lattice
conf=0                    ! counter for conformations

vsol=0.030                ! volume solvent molecule in (nm)^3
vpol= ((4.0/3.0)*pi*(0.3)**3)/vsol  ! volume polymer segment in units of vsol

! eps

eps(1)=eps1
do i=2,ntot
eps(i)=0
enddo


!!!!
! solver
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     init guess all 1.0 

do i=1,n
xg1(i)=1.0
x1(i)=1.0
  do is=1,npoorsv
  xg1(i+n*is)=0.0001
  x1(i+n*is)=0.0001
  enddo
zc(i)= (i-0.5) * delta
enddo

!     init guess from files fort.100 (solvent) and fort.200 (potential)                      

if (infile.ge.1) then
do i=1,n
read(100,*)trash,xfile(i)   ! solvent
x1(i)=xfile(i)
xg1(i)=xfile(i)
  do is=1,npoorsv 
  read(100+is,*)trash,xfile(i+n*is)   ! poorsolvent desde 1 a npoorsv 
  if(xfile(i+n*is).lt.1.0d-30)xfile(i+n*is)=1.0d-30
  x1(i+n*is)=xfile(i+n*is)
  xg1(i+n*is)=xfile(i+n*is)
  enddo !is
enddo !i 
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inn = 0

   call initcha              ! init matrices for chain generation
   conf=0                    ! counter of number of conformations

   do while (conf.lt.cuantas)

   call cadenas(chains,ncha)

   do j=1,ncha

   if(conf.lt.cuantas) then
   conf=conf+1

   do ii = 1, maxntot ! position of first segment


      minpos(conf,ii) = ntot
      maxpos(conf,ii) = 0 

      in1tmp = 0


      do k=1,long
      select case (abs(curvature))
      case (2)
        tempr=((chains(1,k,j)+(float(ii)-0.5)*delta)**2 + chains(2,k,j)**2 +chains(3,k,j)**2 )**(0.5)
        temp=int(tempr/delta)+1  ! put them into the correct layer
      case (1)
        tempr=((chains(1,k,j)+(float(ii)-0.5)*delta)**2+chains(2,k,j)**2)**(0.5)
        temp=int(tempr/delta)+1  ! put them into the correct layer
      case (0) 
        tempr=abs(chains(1,k,j)+(float(ii)-0.5)*delta)
        temp=int(tempr/delta)+1  ! put them into the correct layer
        
      endselect
        
       if(temp.gt.ntot) then
        if(rank.eq.0)print*, 'Increase ntot'
        stop
       endif


       in1tmp(k) = temp

       if(temp.lt.minpos(conf,ii))minpos(conf,ii)=temp
       if(temp.gt.maxpos(conf,ii))maxpos(conf,ii)=temp
       enddo ! k

       if((maxpos(conf,ii)-minpos(conf,ii)).ge.base) then
       print*,'Rank', rank, 'Increase base'
       call MPI_FINALIZE(ierr) ! finaliza MPI
       stop
       endif

       do k = 1, long
       temp = in1tmp(k)-minpos(conf,ii)+1 
       inn(segpoorsv(k),conf,ii,temp) = inn(segpoorsv(k),conf,ii,temp) + 1      
       enddo

   enddo ! ii
   endif

   enddo ! j
   enddo ! while

if(rank.eq.0)print*," chains ready"

! CHECK that chains are unbiased

!   av1 = 0
!   av2 = 0

!   do ii = 1, maxntot
!   do i = 1, cuantas
!   do k = 1, ntot
!   av1(k) = av1(k) + in1n(i,ii,k)
!   av2(k) = av2(k) + in2n(i,ii,k)
!   enddo
!   enddo

!   do k = 1, ntot
!   print*, ii, k, av1(k), av2(k)
!   enddo
!   enddo
!   stop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     computation starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     initializations of input depended variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


iter=0                    ! iteration counter

!if(rank.eq.0) then
open(unit=533,file='lnq.dat')
!open(unit=534,file='ADS-cad.nm-2.dat')
!open(unit=535,file='meanz.dat')
!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MAIN LOOP OVER LAYERS 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



countfile=1


actionflag = 0 ! Actionflag controls the current action of loop
               ! = 0 loop maxntotcounter from 2 to maxntot
               ! = 1 increases npol from npolini to npollast
               ! = 2 decreases npol from npolini to npolfirst
               ! = 3 finalize
maxntotcounter = maxntotcounter_ini !maxntot inicial

npol = npolini

do while (actionflag.lt.3)

 123 if(rank.eq.0)print*, ' npol:', npol, 'maxntot:', maxntotcounter


! xh bulk
 xsolbulk=1.0

do i=1,(npoorsv+1)*n             ! initial guess for x1
xg1(i)=x1(i)
enddo


do i=1,n
do is=1,npoorsv
  if(xg1(i+n*is).lt.1.0d-30)xg1(i+n*is)=1.0d-30 ! OJO
  enddo
enddo


! JEFE
if(rank.eq.0) then ! solo el jefe llama al solver
   iter = 0
   print*, 'solve: Enter solver ', (npoorsv+1)*ntot, ' eqs'
   call call_kinsol(x1, xg1, ier)
   flagsolver = 0
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
endif
! Subordinados


if(rank.ne.0) then
  do
     flagsolver = 0
     source = 0
     CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
     if(flagsolver.eq.1) then
        call call_fkfun(x1) ! todavia no hay solucion => fkfun 
     endif ! flagsolver
     if(flagsolver.eq.0) exit ! Detiene el programa para este nodo
   enddo
endif


! Recupero el valor de ier y de la norma
! Asi los subordinados se enteran si el solver convergio o si hay que
! cambiar la   estrategia...
! Jefe

if (rank.eq.0) then
   norma_tosend = norma
   ier_tosend = ier ! distinto tipo de integer
   CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
   CALL MPI_BCAST(ier_tosend,1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
endif

! Subordinados

if (rank.ne.0) then
   CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
   CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
   norma = norma_tosend
   ier = ier_tosend
endif


do i=1,n
xsol(i)=x1(i)
enddo


if((norma.gt.error).or.(ier.lt.0).or.(isnan(norma))) then
  if(actionflag.gt.0)stop
!stop
!if(rank.eq.0)print*, 'Fail', npol
!if(ccc.eq.1) then
!npol = npol/2.0
!if(rank.eq.0)print*, 'Try', npol
!x1 = xg1
!goto 123
!endif
!npol=(npols(ccc-1)+npol)/2.0
!if(rank.eq.0)print*, 'Try', npol
!x1 = xg1
!goto 123
endif

print*, "norma", rank

!if(npols(ccc).ne.npol) then
!npols(ccc-1) = npol
!npol = npols(ccc)
!goto 123
!endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calc free energy
!
!

call calc_free_energy(actionflag, countfile)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Determination of adsorbed polymer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0) then


write(533,*) npol, dlog(xpol(1))-dlog(q(1))
flush(533)

write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', actionflag,'.',countfile,'.dat'

do is=0,Npoorsv
write(denspolfilename(is),'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densitypolymer',is,'.',actionflag,'.',countfile,'.dat'
enddo

write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', actionflag,'.',countfile,'.dat'
write(totalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densitytotal.',actionflag,'.',countfile,'.dat'
write(xtotalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'xdensitytota.',actionflag,'.',countfile,'.dat'

write(lnqfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'chemical_potent.',actionflag,'.',countfile,'.dat'

open(unit=310,file=sysfilename)
do is=0,Npoorsv
open(unit=1320+is,file=denspolfilename(is))
enddo
open(unit=328,file=totalfilename)
open(unit=329,file=xtotalfilename)
open(unit=330,file=denssolfilename)

open(unit=324,file=lnqfilename)


avtmp=0
do i=1,n
  do is=0,Npoorsv
  write(1320+is,*)zc(i),avpol(is,i)
  avtmp = avtmp + avpol(is,i)
  enddo
write(328,*)zc(i),avtmp
write(329,*)zc(i),xpol(i)
write(330,*)zc(i),xsol(i)
enddo

do i = 1, maxntot
write(324,*)zc(i),dlog(xpol(i))-dlog(q(i))
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     additional system information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(310,*)'GIT Version: ', _VERSION
write(310,*)'system      = neutral polymer'
write(310,*)'fnorm       = ', norma ! residual size of iteration vector
write(310,*)'error       = ',error
write(310,*)'q           = ',q
write(310,*)'length seg  = ',0.50 ! value see subroutine cadenas
write(310,*)'delta       = ',delta
write(310,*)'vsol        = ',vsol
write(310,*)'vpol        = ',vpol*vsol

write(310,*)'npol       = ', npol
write(310,*)'st         = ', st
write(310,*)'Actionflag = ', actionflag

write(310,*)'cuantas     = ',cuantas
write(310,*)'iterations  = ',iter


close(310)
CLOSE(324)
do is=0,Npoorsv
CLOSE(1320+is)
enddo
CLOSE(328)
CLOSE(329)
close(330)

print*, rank, " escribe"

endif ! rank


countfile = countfile+1 ! next
!npolini, npolfirst, npollast, npolstep
select case (actionflag)
 case(0) ! maxntot loop
!    countfile = 1
    if (maxntotcounter.lt.maxntot) then
!    write(1000+maxntotcounter,*), maxntotcounter
    maxntotcounter=maxntotcounter+1
!    xg1 = x1
    endif
    if(maxntotcounter.eq.maxntot) then
    actionflag=1 
    countfile = 1
    endif
 case(1)  ! increases from npolini to npollast
    if(npol.eq.npolini)x1ini=x1
    npol = npol + npolstep      
    if(npol.gt.npollast) then
       npol = npolini
       actionflag = 2
       countfile = 1
       x1 = x1ini
    endif       
 case(2)
    npol = npol - npolstep
    if(npol.lt.npolfirst)actionflag = 3
endselect

END do ! loop de npol

close(533)

call MPI_FINALIZE(ierr) ! finaliza MPI
stop

end

