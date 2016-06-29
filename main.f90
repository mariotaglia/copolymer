
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

implicit none
integer *4 ier ! Kinsol error flag
real*8 pi
real*8 Na               
parameter (Na=6.02d23)

real*8 xsol(ntot)         ! volume fraction solvent
real*8 avtmp
real*8 x1(2*ntot),xg1(2*ntot)   ! density solvent iteration vector
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


integer il,inda,ncha

REAL*8 xfile(2*ntot)                        
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
character*27 denspol2filename
character*27 denspol1filename

integer countfile         ! enumerates the outputfiles 
integer countfileuno     ! enumerates the outputfiles para una corrida
integer conf              ! counts number of conformations

integer readsalt          !integer to read salt concentrations


INTEGER cc, ccc

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

integer in1tmp(long)

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
xg1(i+n)=0.0001
x1(i+n)=0.0001
zc(i)= (i-0.5) * delta
enddo

!     init guess from files fort.100 (solvent) and fort.200 (potential)                      

if (infile.ge.1) then
do i=1,n
read(100,*)trash,xfile(i)   ! solvent
read(200,*)trash,xfile(i+n)   ! solvent
x1(i)=xfile(i)
xg1(i)=xfile(i)
if(xfile(i+n).eq.0.0)xfile(i+n)=1.0d-30
x1(i+n)=xfile(i+n)
xg1(i+n)=xfile(i+n)
enddo  


endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

in1n = 0
in2n = 0

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
        tempr=(chains(1,k,j)+(float(ii)-0.5)*delta)
        temp=int(tempr/delta)+1  ! put them into the correct layer
        
        if(temp.le.0)temp=-temp+2      ! RBC for planar calculation

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



       if(k.le.(long/2))in1n(conf,ii,temp) =  in1n(conf,ii,temp) + 1
       if(k.gt.(long/2))in2n(conf,ii,temp) =  in2n(conf,ii,temp) + 1
!        in1n(conf,ii,temp) =  in1n(conf,ii,temp) + 1
       enddo
   enddo ! ii
   endif

   enddo ! j
   enddo ! while

if(rank.eq.0)print*," chains ready"


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


countfileuno=1           

countfile=1

do cc = 1, nst !loop st
st = sts(cc)

do ccc = 1, nnpol !loop kbind
npol = npols(ccc)

if(rank.eq.0)print*, 'st:',st,' npol:', npol

! xh bulk
 123 xsolbulk=1.0

do i=1,2*n             ! initial gues for x1
xg1(i)=x1(i)
enddo

! JEFE
if(rank.eq.0) then ! solo el jefe llama al solver
   iter = 0
   print*, 'solve: Enter solver ', 2*ntot, ' eqs'
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

do i=1,n
xsol(i)=x1(i)
enddo

if(norma.gt.error) then
if(ccc.eq.1) then
npol = npol/2.0
goto 123
endif
if(rank.eq.0)print*, 'Fail', npol
npol=(npols(ccc-1)+npol)/2.0
if(rank.eq.0)print*, 'Try', npol
npols(ccc-1) = npol
goto 123
endif

if(npols(ccc).ne.npol) then
npol = npols(ccc)
goto 123
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determination of adsorbed polymer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0) then


write(533,*)st, npol, -dlog(qall)

write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', countfileuno,'.',countfile,'.dat'
write(denspol1filename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitypolymer1.',countfileuno,'.',countfile,'.dat'
write(denspol2filename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitypolymer2.',countfileuno,'.',countfile,'.dat'
write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', countfileuno,'.',countfile,'.dat'
write(totalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densitytotal.',countfileuno,'.',countfile,'.dat'
write(xtotalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'xdensitytota.',countfileuno,'.',countfile,'.dat'

write(denspol1filename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitypolymer1.',countfileuno,'.',countfile,'.dat'
write(lnqfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'chemical_potent.',countfileuno,'.',countfile,'.dat'

open(unit=310,file=sysfilename)
open(unit=321,file=denspol1filename)
open(unit=322,file=denspol2filename)
open(unit=323,file=totalfilename)
open(unit=325,file=xtotalfilename)
open(unit=330,file=denssolfilename)
open(unit=324,file=lnqfilename)

do i=1,n
write(321,*)zc(i),avpol(i,1)
write(322,*)zc(i),avpol(i,2)
avtmp = avpol(i,2)+avpol(i,1)
write(323,*)zc(i),avtmp
write(325,*)zc(i),xpol(i)
write(324,*)zc(i),dlog(xpol(i))-dlog(q(i))
write(330,*)zc(i),xsol(i)
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


write(310,*)'cuantas     = ',cuantas
write(310,*)'iterations  = ',iter


close(310)
close(320)
CLOSE(321)
CLOSE(322)
CLOSE(323)
CLOSE(325)
CLOSE(324)
close(330)

countfile = countfile+1 ! next

endif ! rank

END do ! loop de npol
end do ! loop de st

close(533)

countfileuno = countfileuno + 1

call MPI_FINALIZE(ierr) ! finaliza MPI
stop

end

