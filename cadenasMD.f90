!*************************************************************
subroutine cadenasMD(chains,ncha,Uconf, Ntconf,Ugyr, Rgyr, NC)
use seed1
use pis
use matrices
use senos
use globals
use mkai
use longs
use cadenaMD
use MPI
implicit none
integer i,state,j,k1,k2,ncha, is, jj
real*8 rn,dista
real*8 rands,angle
real*8 m(3,3), mm(3,3), m_branch(3,3,50)
real*8 x(3),xend(3,maxlong+5),xendr(3,maxlong+5), xendcom(3,maxlong+5), xend_branch(3,50)
REAL*8 chains(3,maxlong,ncha_max), Uconf
REAL*8 tolerancia    !tolerancia en el calculo de selfavoiding
integer*1 Ntconf(maxlong), seglength(0:Npoorsv)
real*8 Ugyr, Rgyr(0:Npoorsv+1)
real*8 distance(maxlong,maxlong)
integer state_branch(50)
real*8 xendt(3)
integer NC
logical itsopen 
character(80) :: line
integer natoms
integer MDid, MDmol, MDtype
real*8 MDx, MDy, MDz
integer isH
integer timestep


Uconf = 0.0
Ntconf = 0
Ugyr = 0.0
Rgyr = 0.0

! read one conformation
!ITEM: TIMESTEP
!300000
!ITEM: NUMBER OF ATOMS
!93
!ITEM: BOX BOUNDS pp pp pp
!-1.5000000000000000e+02 1.5000000000000000e+02
!-1.5000000000000000e+02 1.5000000000000000e+02
!-1.5000000000000000e+02 1.5000000000000000e+02
!ITEM: ATOMS id mol type x y z


read (7777, '(A)') line
read (7777, *) timestep
read (7777, '(A)') line

print*,'Rank', rank, 'TSTEP', timestep, 'LINEPOS', lineposMD

read (7777,*) natoms
!print*, natoms

do i = 1, 5
read (7777, '(A)') line
!print*,line
enddo

j = 0 ! counts non-H atoms
do i = 1, natoms

  read(7777,*) MDid, MDmol, MDtype, MDx,MDy,MDz

  isH = 0 
  do jj = 1, nMDH
   if(MDHs(jj).eq.MDtype)isH = 1
  enddo

  if(isH.eq.0) then ! store atom
  j = j + 1
   xend(1,j) = MDx/10.0
   xend(2,j) = MDy/10.0
   xend(3,j) = MDz/10.0
  endif 
enddo

if(j.ne.long(NC)) then
  if(rank.eq.0)print*,'Error in MD input file, check long in DEFINITIONS', j
  stop
endif


ncha=0

do i=1,12

  call com(xend,xendcom,long(NC))       ! substracts center of mass
  call rota(xendcom,xendr,long(NC))   ! rotate chain conformation ncha time
  ncha=ncha+1

  if(entflag.eq.1)call print_ent2(xendr,ncha,NC)

  do j=1,long(NC)
    chains(1,j,ncha)=xendr(1,j)       ! output 
    chains(2,j,ncha)=xendr(2,j)
    chains(3,j,ncha)=xendr(3,j)
  enddo
enddo

if(entflag.eq.1)stop

lineposMD = lineposMD + 9 + natoms 

return
end
