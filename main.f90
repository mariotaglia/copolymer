
call initmpi
!call read
call parser
call allocation
call kai
call solve

end

subroutine solve

use globals
use mcharge 
use partfunc
use layer
use volume
use bulk
use seed1
use longs
use MPI
use mkai
use transgauche

implicit none

real*8 solvetime1, solvetime2, solveduration
integer is,ic
integer *4 ier ! Kinsol error flag
real*8 pi
real*8 Na               
parameter (Na=6.02d23)
integer av1(ntot), av2(ntot)
real*8 avtmp
real*8 x1((npoorsv+2)*ntot),xg1((npoorsv+2)*ntot),x1ini((npoorsv+2)*ntot)   ! density solvent iteration vector
real*8 zc(ntot)           ! z-coordinate layer 

REAL*8 sumrhoz, meanz     ! Espesor medio pesado
real*8 pro                ! probability distribution function 
real*8 trash

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 fnorm              ! L2 norm of residual vector function fcn

external fcnelect         ! function containing the SCMFT eqs for solver
integer i,j,k,m,ii,flag,c, jj, iR,iZ ! dummy indice0s

INTEGER temp_R, temp_Z
real*8 tempr_R, tempr_Z
real*8 tmp

real*8 min1               ! variable to determine minimal position of chain
integer qqq,www,eee

integer il,inda,ncha

REAL*8 xfile((npoorsv+2)*ntot)                        
real*8 algo, algo2                  

real*8 chains(3,long,ncha_max) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
real*8 zp(long)
real*8 Uconf
integer*1 Ntconf(long)
real*8 sum,sumel          ! auxiliary variable used in free energy computation  
real*8 sumpi,sumrho,sumrhopol, sumrho2, sumrho2mol !suma de la fraccion de polimero
real*8 sumUgyr, sumRgyr(0:Npoorsv+1), Rgyr(0:Npoorsv+1), Ugyr, Rgyrprom(0:Npoorsv+1)

! global running files
!character*15 meanzfilename
!character*15 sigmafilename
!character*17 sigmaadfilename

! single layer files
character*18 sysfilename      ! contains value of free energy, input parameter etc
character*29 phifilename      ! electric potential 
character*26 denssolfilename  ! contains the denisty of the solvent
character*27 lnqfilename  ! contains the denisty of the solvent
character*28 densendfilename
CHARACTER*24 totalfilename
CHARACTER*24 xtotalfilename
CHARACTER*18 ntransfilename
character*27 densposfilename
character*27 dielfilename
character*27 densnegfilename
character*24 densHplusfilename
character*24 densOHminfilename
character*50, allocatable :: denspolfilename(:)
character*50, allocatable :: fracBHplus(:)
character*48, allocatable :: fracAmin(:)
character*47, allocatable :: densAcidfilename(:)
character*48, allocatable :: densBasicfilename(:)

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

integer in1tmp(long,2)

allocate(denspolfilename(0:Npoorsv))
allocate(fracAmin(Nacids),fracBHplus(Nbasics))
allocate(densAcidfilename(Nacids),densBasicfilename(Nbasics))

error = 1.0d-6
!seed=435 !+ 3232*rank               ! seed for random number generator

if(rank.eq.0) then
  seed=435
  print*, 'I am', rank, ' and my seed is', seed
else 
  seed=0
endif

if(rank.eq.0)print*, 'Program Multicapa'
if(rank.eq.0)print*, 'GIT Version: ', _VERSION

!     initializations of variables 
pi=dacos(-1.0d0)          ! pi = arccos(-1) 
itmax=200                 ! maximum number of iterations       
n=ntot                    ! size of lattice
conf=0                    ! counter for conformations

vsol=0.030                ! volume solvent molecule in (nm)^3
vpol(:)=vpol(:)/vsol  ! volume polymer segment in units of vsol
vpol_a(:)=vpol_a(:)/vsol
vpol_b(:)=vpol_b(:)/vsol

vchain=0.0
do i=1,long
  vchain=vchain+vpol(segpoorsv(i))
enddo

vneg=4/3*pi*0.2**3/vsol !volume of anion in units of vsol
vpos=4/3*pi*0.2**3/vsol !volume of cation in units of vsol 

pKw=14.0

cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
rhosalt=Csalt*Na/(1.0d24) !salt conc. in unit of nº of particles/nm³
wperm = 0.114 !water permitivity in units of e^2/kT.nm

if(pHbulk.le.7) then  ! pH<= 7
  xposbulk= rhosalt*vsol*vpos
  xnegbulk= rhosalt*vsol*vneg + (xHplusbulk-xOHminbulk)*vneg ! NaCl+ HCl  
else                  ! pH >7 
  xposbulk= rhosalt*vsol*vpos + (xOHminbulk-xHplusbulk)*vpos ! NaCl+ NaOH   
  xnegbulk= rhosalt*vsol*vneg
endif


xsolbulk=1-xposbulk-xnegbulk-xHplusbulk-xOHminbulk ! bulk volume fraction of solvent 

if (Nacids.gt.0) then
  do i=1,Nacids
    Ka(i) = 10**(-pKa(i))
    Ka(i) = (Ka(i)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant
  enddo
endif

if (Nbasics.gt.0) then
  do i=1,Nbasics
    Kb(i) = 10**(-pKb(i))
    Kb(i) = (Kb(i)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant
  enddo
endif

!expmupos=rhosalt*vsol*vpos/xsolbulk**vpos
expmupos=xposbulk/vpos/xsolbulk**vpos  
!expmuneg=rhosalt*vsol*vneg/xsolbulk**vneg 
expmuneg=xnegbulk/vneg/xsolbulk**vneg           

expmuHplus=xHplusbulk/xsolbulk ! vHplus=vsol
expmuOHmin=xOHminbulk/xsolbulk ! vOHminus=vsol


print*, "I am rank", rank, "and I generate and calculate", cuantas, "conformation out of", totalcuantas

!!!!
! solver
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     init guess all 1.0 

do iR=1,dimR
do iZ=1,dimZ

   xg1(dimR*(iZ-1)+iR)=1.0
   x1(dimR*(iZ-1)+iR)=1.0

   do is=1,Npoorsv+1
      xg1(n*is+dimR*(iZ-1)+iR)=0.0001
      x1(n*is+dimR*(iZ-1)+iR)=0.0001
   enddo

!   zc(iZ)= (iZ-0.5) * delta !OJO
enddo
   zc(iR)= (iR-0.5) * deltaR
enddo

!     init guess from files fort.100 (solvent) and fort.200 (potential)                      

if (infile.ge.1) then

   do iR=1,dimR
   do iZ=1,dimZ

     read(100,*)trash,trash,xfile(dimR*(iZ-1)+iR)   ! solvent
     x1(dimR*(iZ-1)+iR)=xfile(dimR*(iZ-1)+iR)
     xg1(dimR*(iZ-1)+iR)=xfile(dimR*(iZ-1)+iR)

     do is=1,npoorsv 

       read(100+is,*)trash,trash,xfile(n*is+dimR*(iZ-1)+iR)   ! poorsolvent desde 1 a npoorsv 
!       if(xfile(i+n*is).lt.1.0d-30)xfile(i+n*is)=1.0d-30
       x1(n*is+dimR*(iZ-1)+iR)=xfile(n*is+dimR*(iZ-1)+iR)
       xg1(n*is+dimR*(iZ-1)+iR)=xfile(n*is+dimR*(iZ-1)+iR)

     enddo !is

     read(200,*)trash,trash,xfile(n*(Npoorsv+1)+dimR*(iZ-1)+iR)
     x1(n*(npoorsv+1)+dimR*(iZ-1)+iR)=xfile(n*(npoorsv+1)+dimR*(iZ-1)+iR)
     xg1(n*(npoorsv+1)+dimR*(iZ-1)+iR)=xfile(n*(npoorsv+1)+dimR*(iZ-1)+iR)

   enddo !iR
   enddo !iZ

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call initcha              ! init matrices for chain generation

conf=0                    ! counter of number of conformations

if (rank.gt.0) then
   call MPI_RECV(seed, 1, MPI_INTEGER, rank-1, rank-1, MPI_COMM_WORLD, status, ierr)
   print*, "rank", rank, "received from rank", rank-1,"seed =", seed
endif

innR = 0
innZ = 0
sumRgyr(:)=0.
sumUgyr=0.
Rgyrprom(:)=0.

do while (conf.lt.cuantas)

!   print*, "conf is",conf,". Seed is", seed
  
   call cadenas(chains,ncha,Uconf,Ntconf,Ugyr,Rgyr)
   
   do is=0,Npoorsv+1
      sumRgyr(is)=sumRgyr(is)+Rgyr(is)*exp(-Ugyr)
   enddo

   sumUgyr=sumUgyr+exp(-Ugyr)

   do j=1,ncha

      if(conf.lt.cuantas) then

         conf=conf+1
         Uchain(conf)=Uconf
         Ntrans(:,conf) = Ntconf(:)

         do k=1,long
            do ii = 1,maxntotR ! position of first segment (or Center of mass?)

               select case (abs(curvature))
                 case (2)
                  tempr_R=((chains(1,k,j)+(float(ii)-0.5)*deltaR)**2 + chains(2,k,j)**2 + chains(3,k,j)**2 )**(0.5)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case (1)
                  tempr_R=((chains(1,k,j)+(float(ii)-0.5)*deltaR)**2+chains(2,k,j)**2)**(0.5)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case (0) 
                  tempr_R=abs(chains(1,k,j)+(float(ii)-0.5)*deltaR)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
               endselect
             
               if(temp_R.gt.dimR) then
                  if(rank.eq.0)print*, 'main.f90: increase dimR'
                  stop
               endif

              innR(k,conf,ii)=temp_R ! in which layer is the segment "k" of a chain at position "ii" and conformation "conf"
         
            enddo ! ii
         
         tempr_Z=chains(3,k,j)
         innZ(k,conf)=int(anint(tempr_Z/deltaZ))

         enddo ! k

      endif

   enddo ! j

enddo ! while

if (rank.lt.size-1) then
   call MPI_SEND(seed, 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD, ierr)
   print*, "rank", rank, "sent to rank", rank+1,"seed=", seed
endif

do is=0,Npoorsv+1
   Rgyrprom(is)=sumRgyr(is)/sumUgyr
enddo

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if(rank.eq.0) then

   print*," chains ready"

   do is=0,Npoorsv+1
      print*,is, Rgyrprom(is), sumRgyr(is), sumUgyr
   enddo

   !do k = 1, 20
      !print*,10*k,Uchain(10*k)
   !enddo

endif

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
open(unit=534,file='rgyration.dat')
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

! maxntotcounter = maxntotcounter_ini !maxntot inicial
maxntotcounterR = maxntotR
maxntotcounterZ = maxntotZ
npol = npolini

do while (actionflag.lt.3)

   123 if(rank.eq.0)print*, ' npol:', npol, 'maxntotR:', maxntotcounterR, 'maxntotZ:', maxntotcounterZ


! xh bulk
! xsolbulk=1.0

   do i=1,(npoorsv+2)*n             ! initial guess for x1
      xg1(i)=x1(i)
   enddo

!   do i=1,n
!     do is=1,npoorsv+1
!        if(xg1(i+n*is).lt.1.0d-30)xg1(i+n*is)=1.0d-30 ! OJO
!     enddo
!   enddo

! JEFE
   if(rank.eq.0) then ! solo el jefe llama al solver
      print*, 'solve: Enter solver ', (npoorsv+2)*ntot, ' eqs'
   endif   

   iter=0
   call call_kinsol(x1, xg1, ier)
!      flagsolver = 0
!   endif

!   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)

!   print*, "I am rank", rank, "and I receive flagsolver"

!   call MPI_BARRIER
   
!   Subordinados

!   if(rank.ne.0) then

!      do

!         flagsolver = 1
!         source = 0

!         if(flagsolver.eq.1) then
!           call call_fkfun(x1) ! todavia no hay solucion => fkfun 
!         endif ! flagsolver

!         if(flagsolver.eq.0) exit ! Detiene el programa para este nodo

!      enddo

!   endif


! Recupero el valor de ier y de la norma
! Asi los subordinados se enteran si el solver convergio o si hay que
! cambiar la   estrategia...
! Jefe

!   if (rank.eq.0) then
!      norma_tosend = norma
!      ier_tosend = ier ! distinto tipo de integer
!   endif
!      CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!      CALL MPI_BCAST(ier_tosend,1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
!   endif

! Subordinados

!   if (rank.ne.0) then
!      CALL MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
!      CALL MPI_BCAST(ier_tosend, 1, MPI_INTEGER,0,MPI_COMM_WORLD,err)
!      norma = norma_tosend
!      ier = ier_tosend
!   endif


   do iR=1,dimR
   do iZ=1,dimZ
     xsol(iR,iZ)=x1(dimR*(iZ-1)+iR)
   enddo
   enddo

   if((norma.gt.error).or.(ier.lt.0).or.(isnan(norma))) then
      if(actionflag.gt.0) then
         print*, " I am ", rank, " I stopped the work"
         stop
      endif
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

      write(533,*) npol, dlog(xpol(1,1))-dlog(q(1,1))
      flush(533)

      write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', actionflag,'.',countfile,'.dat'

      write(phifilename,'(A18,BZ,I3.3,A1,I3.3,A4)')'electricpotential.', actionflag,'.',countfile,'.dat'

      do is=0,Npoorsv
         write(denspolfilename(is),'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densitypolymer',is,'.',actionflag,'.',countfile,'.dat'
      enddo

      if (Nacids.ge.1) then
        do ic=1,Nacids
          write(fracAmin(ic),'(A12,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'fractionAmin',ic,'.',actionflag,'.',countfile,'.dat'
          write(densAcidfilename(ic),'(A11,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densityacid',ic,'.',actionflag,'.',countfile,'.dat'
        enddo
      endif

      if (Nbasics.ge.1) then
        do ic=1,Nbasics
          write(fracBHplus(ic),'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'fractionBHplus',ic,'.',actionflag,'.',countfile,'.dat'
          write(densBasicfilename(ic),'(A12,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densitybasic',ic,'.',actionflag,'.',countfile,'.dat'
        enddo
      endif

      write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', actionflag,'.',countfile,'.dat'
      write(densposfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitypositive.',actionflag,'.',countfile,'.dat'
      write(dielfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'dielectric_cons.',actionflag,'.',countfile,'.dat'
      write(densnegfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitynegative.',actionflag,'.',countfile,'.dat'
      write(densHplusfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densityHplus.',actionflag,'.',countfile,'.dat'
      write(densOHminfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densityOHmin.',actionflag,'.',countfile,'.dat'

      write(totalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densitytotal.',actionflag,'.',countfile,'.dat'

      write(xtotalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'xdensitytota.',actionflag,'.',countfile,'.dat'
      write(ntransfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'ntrans.',actionflag,'.',countfile,'.dat'

      write(lnqfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'chemical_potent.',actionflag,'.',countfile,'.dat'

      open(unit=310,file=sysfilename)
      open(unit=311,file=phifilename)

      do is=0,Npoorsv
         open(unit=1320+is,file=denspolfilename(is))
      enddo
       
      do ic=1,Nacids
        open(unit=1050+ic, file=fracAmin(ic))
        open(unit=1780+ic, file=densAcidfilename(ic))
      enddo
  
      do ic=1,Nbasics
        open(unit=1520+ic, file=fracBHplus(ic))
        open(unit=1680+ic, file=densBasicfilename(ic))
      enddo

      open(unit=328,file=totalfilename)
      open(unit=329,file=xtotalfilename)
      open(unit=327,file=ntransfilename)
      open(unit=330,file=denssolfilename)
      open(unit=331,file=densposfilename)
      open(unit=332,file=densnegfilename)
      open(unit=333,file=densHplusfilename)
      open(unit=334,file=densOHminfilename)
      open(unit=324,file=lnqfilename)
      open(unit=335,file=dielfilename)

      do i = 3, long-1
         write(327,*)i, trans(i)
      enddo

      avtmp=0

      do iR=1,dimR
      do iZ=1,dimZ

         do is=0,Npoorsv
            write(1320+is,*)zc(iR),iZ,avpol(is,iR,iZ)
            avtmp = avtmp + avpol(is,iR,iZ)
         enddo

         if (Nacids.ge.1) then
           do ic=1,Nacids
             write(1050+ic,*)zc(iR),iZ,fAmin(ic,iR,iZ)
             write(1780+ic,*)zc(iR),iz,avpola(ic,iR,iZ)
           enddo
         endif
       
         if (Nbasics.ge.1) then
           do ic=1,Nbasics
             write(1520+ic,*)zc(iR),iZ,fBHplus(ic,iR,iZ)
             write(1680+ic,*)zc(iR),iZ,avpolb(ic,iR,iZ)
           enddo
         endif
 
         write(328,*)zc(iR),iZ,avtmp
         write(329,*)zc(iR),iZ,xpol(iR,iZ)
         write(330,*)zc(iR),iZ,xsol(iR,iZ)
         write(331,*)zc(iR),iZ,avpos(iR,iZ)
         write(332,*)zc(iR),iZ,avneg(iR,iZ)
         write(333,*)zc(iR),iZ,avHplus(iR,iZ)
         write(334,*)zc(iR),iZ,avOHmin(iR,iZ)
         write(335,*)zc(iR),iZ,epsfcn(iR,iZ)
         write(311,*)zc(iR),iZ,phi(iR,iZ)

      enddo
      enddo

      do iR = 1, maxntotR
      do iZ = 1, maxntotZ
          write(324,*)zc(i),iR,dlog(xpol(iR,iZ))-dlog(q(iR,iZ))
      enddo
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
      write(310,*)'deltaR      = ',deltaR
      write(310,*)'deltaZ      = ',deltaZ
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

      if (Nacids.ge.1) then  
        do ic=1,Nacids
          close(1050+ic)
          close(1780+ic)
        enddo
      endif
  
      if (Nbasics.ge.1) then
        do ic=1,Nbasics
          close(1520+ic)
          close(1680+ic)
        enddo
      endif

      CLOSE(327)
      CLOSE(328)
      CLOSE(329)
      close(330)
      close(331)
      close(332)
      close(333)
      close(334)
      close(311)

      print*, rank, " escribe"

   endif ! rank

   countfile = countfile+1 ! next
   !npolini, npolfirst, npollast, npolstep

   select case (actionflag)

     case(0) ! maxntot loop

!     countfile = 1
      if (maxntotcounterR.lt.maxntotR) then
!       write(1000+maxntotcounter,*), maxntotcounter
        maxntotcounterR=maxntotcounterR+1
!       xg1 = x1
      endif

      if(maxntotcounterR.eq.maxntotR) then
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

end do ! loop de npol

do is=0,Npoorsv+1
   write(534,*) is, Rgyrprom(is)
enddo

close(533)
close(534)

!if (rank.eq.0) then
!   solvetime2=MPI_WTIME()
!   solveduration=solvetime2-solvetime1
!   print*, "The solver took", solveduration, "seconds."
!endif

call MPI_FINALIZE(ierr) ! finaliza MPI

stop

end
