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
real*8 tmp
integer counter, counter2
integer is,ic
integer NC
external fcnelect         ! function containing the SCMFT eqs for solver
integer i,iR,iZ ! dummy indice0s
real*8 zc(dimR)
! single layer files
character*50 sysfilename      ! contains value of free energy, input parameter etc
character*50 phifilename      ! electric potential 
character*50 denssolfilename  ! contains the denisty of the solvent
character*50 lnqfilename  ! contains the denisty of the solvent
CHARACTER*50 xtotalfilename
CHARACTER*50 poorsvfilename
CHARACTER*50 ntransfilename
character*50 densposfilename
character*50 denstotfilename
character*50 dielfilename
character*50 densnegfilename
character*50 densHplusfilename
character*50 densOHminfilename
character*50 denspolfilename
character*50 fracBHplus
character*50 fracAmin
character*50 densAcidfilename
character*50 densBasicfilename
character*50 outfilename


do iR = 1, dimR
zc(iR)= (iR-0.5) * deltaR
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiple run files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do NC = 1, Ncomp
write(1533+NC,*) npol, dlog(xpol(1,1,NC))-dlog(q(1,1,NC))
flush(1533+NC)
enddo

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
      write(310,*)'rsalt = ',r_pos,r_neg
      write(310,*)'dielP = ',dielP
      write(310,*)'npol       = ', npol
      write(310,*)'st         = ', st
      write(310,*)'Actionflag = ', counter

      write(310,*)'cuantas     = ',cuantas
      write(310,*)'iterations  = ',iter


      write(310,*)'npolratios  = ',npolratio

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

do NC=1, Ncomp
do is=0,Npoorsv

write(denspolfilename,'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A1,I3.3,A4)')'densitypolymer',NC,'.',is,'.',counter,'.',counter2,'.dat'
open(unit=311,file=denspolfilename)

do iR=1,dimR
   do iZ=1,dimZ
       write(311,*)zc(iR),iZ,avpol(is,iR,iZ,NC)
   enddo
enddo

close(311)
enddo ! is
enddo ! NC

! Poor sv (use for input)

do is=1,Npoorsv

write(poorsvfilename,'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'xpoorsvpolymer',is,'.',counter,'.',counter2,'.dat'
open(unit=311,file=poorsvfilename)

do iR=1,dimR
   do iZ=1,dimZ
       write(311,*)zc(iR),iZ,xtotal(is,iR,iZ)
   enddo
enddo

close(311)
enddo ! is

! Density polymer total

do is=0,Npoorsv
write(denstotfilename,'(A14,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'densitytotalpol',is,'.',counter,'.',counter2,'.dat'
open(unit=311,file=denstotfilename)

do iR=1,dimR
   do iZ=1,dimZ
       tmp = 0
       do NC = 1, Ncomp 
       tmp = tmp + avpol(is,iR,iZ,NC)
       enddo
       write(311,*)zc(iR),iZ,tmp
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
do NC = 1, Ncomp
do ic=1,Nacids
write(densAcidfilename,'(A11,BZ,I3.3,A1,I3.3, A1, I3.3,A1,I3.3,A4)')'densityacid',NC,'.',ic,'.',counter,'.',counter2,'.dat'
open(unit=1780, file=densAcidfilename)
do iR=1,dimR
   do iZ=1,dimZ
       write(1780,*)zc(iR),iz,avpola(ic,iR,iZ,NC)
   enddo
enddo
close(1780)
enddo
enddo ! NC

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
do NC = 1, Ncomp
do ic=1,NBasics
write(densBasicfilename,'(A11,BZ,I3.3,A1,I3.3,A1,I3.3,A1,I3.3,A4)')'densitybasic',NC,'.',ic,'.',counter,'.',counter2,'.dat'
open(unit=1780, file=densBasicfilename)
do iR=1,dimR
   do iZ=1,dimZ
       write(1780,*)zc(iR),iz,avpolb(ic,iR,iZ,NC)
   enddo
enddo
close(1780)
enddo
enddo ! NC
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
open(unit=334,file=densOHminfilename)
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

do NC = 1, Ncomp
write(xtotalfilename,'(A13,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'xdensitytota.',NC,'.',counter,'.',counter2,'.dat'
open(unit=329,file=xtotalfilename)
 do iR=1,dimR
   do iZ=1,dimZ
      write(329,*)zc(iR),iZ,xpol(iR,iZ,NC)
   enddo
 enddo
close(329)
enddo

! Number trans
do NC = 1, NComp
write(ntransfilename,'(A7,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'ntrans.',NC,'.',counter,'.',counter2,'.dat'
open(unit=327,file=ntransfilename)
do i = 3, long(NC)-1
    write(327,*)i, trans(i,NC)
enddo
close(327)
enddo ! NC

! lnq

do NC = 1,Ncomp
write(lnqfilename,'(A16,BZ,I3.3,A1,I3.3,A1,I3.3,A4)')'chemical_potent.',NC,'.',counter,'.',counter2,'.dat'
open(unit=324,file=lnqfilename)
do iR = 1, maxntotR
   do iZ = 1, maxntotZ
       write(324,*)zc(iR),iZ,dlog(xpol(iR,iZ,NC))-dlog(q(iR,iZ,NC))
    enddo
enddo
close(324)
enddo ! NC

write(outfilename,'(A4, I3.3,A1,I3.3,A4)')'out.', counter, '.', counter2, '.dat'
open(unit=8, file=outfilename, form='unformatted')
write(8)xflag
close(8)


end
