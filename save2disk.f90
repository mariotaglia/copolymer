subroutine save2disk(counter,counter2)

use globals
use mcharge 
use partfunc
use layer
use volume
use bulk
use longs
use mkai
use transgauche
use pis
implicit none

integer counter, counter2
integer is,ic
real*8 zc(dimR)           ! z-coordinate layer 

external fcnelect         ! function containing the SCMFT eqs for solver
integer i,iR,iZ ! dummy indice0s

! single layer files
character*18 sysfilename      ! contains value of free energy, input parameter etc
character*29 phifilename      ! electric potential 
character*26 denssolfilename  ! contains the denisty of the solvent
character*27 lnqfilename  ! contains the denisty of the solvent
CHARACTER*24 xtotalfilename
CHARACTER*18 ntransfilename
character*27 densposfilename
character*27 dielfilename
character*27 densnegfilename
character*24 densHplusfilename
character*24 densOHminfilename
character*50 denspolfilename
character*50 fracBHplus
character*48 fracAmin
character*47 densAcidfilename
character*48 densBasicfilename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiple run files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(533,*) npol, dlog(xpol(1,1))-dlog(q(1,1))
flush(533)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Single run files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SYSTEM

write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', counter,'.',counter2,'.dat'
open(unit=310,file=sysfilename)

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
      write(310,*)'Actionflag = ', counter

      write(310,*)'cuantas     = ',cuantas
      write(310,*)'iterations  = ',iter

      close(310)
 
! Electrostatic potential

write(phifilename,'(A18,BZ,I3.3,A1,I3.3,A4)')'electricpotential.', counter,'.',counter2,'.dat'
open(unit=311,file=phifilename)
do iR=1,dimR
   do iZ=1,dimZ
         write(311,*)zc(iR),iZ,phi(iR,iZ)
   enddo
enddo
close(311)

! Density polymer

do is=0,Npoorsv
write(denspolfilename,'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densitypolymer',is,'.',counter,'.',counter2,'.dat'
open(unit=1320,file=denspolfilename)
do iR=1,dimR
   do iZ=1,dimZ
       write(1320,*)zc(iR),iZ,avpol(is,iR,iZ)
   enddo
enddo
close(311)
enddo ! is

! Frac A min
do ic=1,Nacids
write(fracAmin,'(A12,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'fractionAmin',ic,'.',counter,'.',counter2,'.dat'
open(unit=1050, file=fracAmin)
do iR=1,dimR
   do iZ=1,dimZ
       write(1050,*)zc(iR),iZ,fAmin(ic,iR,iZ)
   enddo
enddo
close(1050)
enddo

! Density acid
do ic=1,Nacids
write(densAcidfilename,'(A11,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densityacid',ic,'.',counter,'.',counter2,'.dat'
open(unit=1780, file=densAcidfilename)
do iR=1,dimR
   do iZ=1,dimZ
       write(1780,*)zc(iR),iz,avpola(ic,iR,iZ)
   enddo
enddo
close(1780)
enddo

! Frac B min
do ic=1,NBasics
write(fracBHplus,'(A12,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'fractionBHplus',ic,'.',counter,'.',counter2,'.dat'
open(unit=1050, file=fracBHplus)
do iR=1,dimR
   do iZ=1,dimZ
       write(1050,*)zc(iR),iZ,fBHplus(ic,iR,iZ)
   enddo
enddo
close(1050)
enddo

! Density basic
do ic=1,NBasics
write(densBasicfilename,'(A11,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densitybasic',ic,'.',counter,'.',counter2,'.dat'
open(unit=1780, file=densBasicfilename)
do iR=1,dimR
   do iZ=1,dimZ
       write(1780,*)zc(iR),iz,avpolb(ic,iR,iZ)
   enddo
enddo
close(1780)
enddo

! Solvent

write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', counter,'.',counter2,'.dat'
open(unit=330,file=denssolfilename)
do iR=1,dimR
   do iZ=1,dimZ
       write(330,*)zc(iR),iZ,xsol(iR,iZ)
   enddo
enddo
close(330)

! cations

write(densposfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitypositive.',counter,'.',counter2,'.dat'
open(unit=331,file=densposfilename)
 do iR=1,dimR
   do iZ=1,dimZ
       write(331,*)zc(iR),iZ,avpos(iR,iZ)
   enddo
enddo
close(331)

! anions

write(densnegfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'densitynegative.',counter,'.',counter2,'.dat'
open(unit=332,file=densnegfilename)
 do iR=1,dimR
   do iZ=1,dimZ
       write(332,*)zc(iR),iZ,avneg(iR,iZ)
   enddo
enddo
close(332)

! protons

write(densHplusfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densityHplus.',counter,'.',counter2,'.dat'
open(unit=333,file=densHplusfilename)
 do iR=1,dimR
   do iZ=1,dimZ
       write(333,*)zc(iR),iZ,avHplus(iR,iZ)
   enddo
enddo
close(333)

! OH-

write(densOHminfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'densityOHmin.',counter,'.',counter2,'.dat'
 do iR=1,dimR
   do iZ=1,dimZ
       write(334,*)zc(iR),iZ,avOHmin(iR,iZ)
   enddo
enddo
close(334)

! dielectric constant

write(dielfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'dielectric_cons.',counter,'.',counter2,'.dat'
open(unit=335,file=dielfilename)
 do iR=1,dimR
   do iZ=1,dimZ
       write(335,*)zc(iR),iZ,epsfcn(iR,iZ)
   enddo
enddo
close(335)

! xtotal

write(xtotalfilename,'(A13,BZ,I3.3,A1,I3.3,A4)')'xdensitytota.',counter,'.',counter2,'.dat'
open(unit=329,file=xtotalfilename)
 do iR=1,dimR
   do iZ=1,dimZ
      write(329,*)zc(iR),iZ,xpol(iR,iZ)
   enddo
 enddo
close(329)

! Number trans

write(ntransfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'ntrans.',counter,'.',counter2,'.dat'
open(unit=327,file=ntransfilename)
do i = 3, long-1
         write(327,*)i, trans(i)
enddo
close(327)

! lnq

write(lnqfilename,'(A16,BZ,I3.3,A1,I3.3,A4)')'chemical_potent.',counter,'.',counter2,'.dat'
open(unit=324,file=lnqfilename)
do iR = 1, maxntotR
   do iZ = 1, maxntotZ
       write(324,*)zc(i),iR,dlog(xpol(iR,iZ))-dlog(q(iR,iZ))
    enddo
enddo
close(324)


close(534)
end
