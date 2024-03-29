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
use mkai, only: segpoorsv
use mcharge, only : acidtype, basictype
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
integer MDid(2000), MDmol(2000), MDtype(2000)
real*8 MDx(2000), MDy(2000), MDz(2000)
integer isH
integer*8 timestep
integer k,kk

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

print*,'CADENASMD, rank', rank, 'TSTEP', timestep, 'LINEPOS', lineposMD

read (7777,*) natoms
!print*, natoms

do i = 1, 5
read (7777, '(A)') line
!print*,line
enddo

do i = 1, natoms
  read(7777,*) MDid(i), MDmol(i), MDtype(i), MDx(i),MDy(i),MDz(i)
enddo


jj = 0 ! current position in MOLT bead list 
do j = 1, natoms ! current position in MD atom list
 do i = 1, natoms
   if (MDid(i).eq.j) then ! we found atom j
     if(MDHs(MDtype(i),NC).eq.1) then ! atom j is a heavy atom
       jj = jj + 1 ! advance one MOLT bead list
       xend(1,jj) = MDx(i)/10.0
       xend(2,jj) = MDy(i)/10.0
       xend(3,jj) = MDz(i)/10.0
       segpoorsv(jj,NC) = MDsegpoorsv(MDtype(i),NC)
       acidtype(jj,NC) = MDacidtype(MDtype(i),NC)
       basictype(jj,NC) = MDbasictype(MDtype(i),NC)
     endif
    endif
 enddo ! i
enddo ! j

if(jj.ne.long(NC)) then
  if(rank.eq.0)print*,'Error in MD input file, check long in DEFINITIONS', j
  stop
endif

! DEBUG
!do i = 1, long(NC)
! print*, 'rank!', rank,i, segpoorsv(i,NC), acidtype(i,NC), basictype(i,NC)
!enddo



ncha=0

do i=1,nrot(NC)

  call com(xend,xendcom,long(NC))       ! substracts center of mass
  call rota(xendcom,xendr,long(NC))   ! rotate chain conformation ncha time
  
  if (flagreflex(NC).eq.1) then

     ncha=ncha+1

     if(entflag.eq.1)call print_ent2(xendr,ncha,NC)

     do j=1,long(NC)
       chains(1,j,ncha)=xendr(1,j)
       chains(2,j,ncha)=xendr(2,j)
       chains(3,j,ncha)=xendr(3,j)
     enddo

  elseif (flagreflex(NC).eq.2) then
    do k=1,2  ! to reflect the chains in the R coordinate
 
      kk=2*(k-1)-1 ! kk = {-1,1} for k = {1,2} to reflect the chains in the R coordinate

      ncha=ncha+1

      if(entflag.eq.1)call print_ent2(xendr,ncha,NC)

      do j=1,long(NC)

       chains(1,j,ncha)=kk*xendr(1,j)   !kk*xendr(1,j)       ! to reflect the chains in the R coordinate
       chains(2,j,ncha)=xendr(2,j)
       chains(3,j,ncha)=xendr(3,j)

      enddo
     enddo

  else
    print*,"flagreflex must be either 1 or 2"
  stop

  endif

enddo

if(entflag.eq.1)stop

lineposMD = lineposMD + 9 + natoms 

return
end
