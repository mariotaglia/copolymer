subroutine mrrrr(a,b,c)
real*8 a(3,3),b(3,3),c(3,3)

do i=1,3
 do j=1,3
  c(i,j)=0
 enddo
enddo

do i=1,3
 do j=1,3
  do k=1,3
    c(i,j)=c(i,j)+a(i,k)*b(k,j)
  enddo
 enddo
enddo

return
end 

subroutine initcha
use pis
use matrices
use senos
implicit none

pi=acos(-1.0000000e0)
sitheta=sin(68.0*pi/180.0)
cotheta=cos(68.0*pi/180.0)
siphip=sin(120.0*pi/180.0)
cophip=cos(120.0*pi/180.0)

tt(1,1)=cotheta
tt(1,2)=sitheta
tt(1,3)=0.0
tt(2,1)=sitheta
tt(2,2)=-cotheta
tt(2,3)=0.0
tt(3,1)=0.0
tt(3,2)=0.0
tt(3,3)=-1.0

tp(1,1)=cotheta
tp(1,2)=sitheta
tp(1,3)=0.0
tp(2,1)=sitheta*cophip
tp(2,2)=-cotheta*cophip
tp(2,3)=siphip
tp(3,1)=sitheta*siphip
tp(3,2)=-cotheta*siphip
tp(3,3)=-cophip

tm(1,1)=cotheta
tm(1,2)=sitheta
tm(1,3)=0.0
tm(2,1)=sitheta*cophip
tm(2,2)=-cotheta*cophip
tm(2,3)=-siphip
tm(3,1)=-sitheta*siphip
tm(3,2)=cotheta*siphip
tm(3,3)=-cophip
return
end


!**************************************************************
!* rotates a given chains conformation                        *  
!* pre: xend = input chain                                    *
!* post: xendr= rotated chain                                 *
!**************************************************************     

subroutine rota(xend,xendr,n)
use seed1
use layer
use pis
implicit none

integer n
real*8 xend(3,n+5),rands,xendr(3,n+5)
integer i
real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
real*8 alfa,gama,cga,a,b,c

fac=rands(seed)
fac1=rands(seed)
fac2=rands(seed)
alfa=fac*2*pi
cbe=fac1*2.0-1.0
gama=fac2*2*pi

sbe=(1-cbe**2)**0.5
cal=cos(alfa)
sal=sin(alfa)
cga=cos(gama)
sga=sin(gama)

do i=1,n+1               ! rotation segmentos

  a=xend(1,i)
  b=xend(2,i)
  c=xend(3,i)

  xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
  xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
  xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe

enddo

return

end


!*************************************************************
!* Generates a chains conformations bases                    *
!* on a three state RIS-model see Flory book                 *
!* GENERA CADENAS DE PAH-Os                                  *
!*************************************************************
subroutine cadenas(chains,ncha,Uconf, Ntconf,Ugyr, Rgyr, NC)
use seed1
use pis
use matrices
use senos
use globals
use mkai
use longs
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

tolerancia = 1.0e-5


223 Uconf=0.0
Ntconf(:) = 0

seglength=0
Ugyr=0.0
Rgyr=0.0
distance(:,:)=0.0
 
xend(1,1)=0.0      ! first position 
xend(2,1)=0.0
xend(3,1)=0.0
rn=rands(seed)
angle=0.0

! transition probalibity
m(1,1)=cotheta          
m(1,2)=sitheta
m(1,3)=0.0
m(2,1)=cos(angle)*sitheta
m(2,2)=-cos(angle)*cotheta
m(2,3)=sin(angle)
m(3,1)=sin(angle)*sitheta
m(3,2)=-sin(angle)*cotheta
m(3,3)=-cos(angle)

x(1)=m(1,1)*lseg     
x(2)=m(2,1)*lseg
x(3)=m(3,1)*lseg

xend(1,2)=xend(1,1)+x(1)  ! second postion
xend(2,2)=xend(2,1)+x(2)
xend(3,2)=xend(3,1)+x(3)


do i=3,long(NC)-long_branches(NC)          ! loop over remaining positions!

rn=rands(seed)
state=int(rn*3)        ! random select the state= {trans,gauch+,gauch-}

  if (state.eq.3) then 
    state=2
  endif

  if (state.eq.0) then
!*********************************** TRANS     
    call mrrrr(m,tt,mm)
    if(i.gt.3)Uconf=Uconf+Ut(segpoorsv(i-1,NC)) ! first segment to have a dihedral angle is i = 4
                                           ! OJO : the order of chain grown changes the assigment of diehdral angles
    if(i.gt.3)Ntconf(i-1) = 1   

  elseif (state.eq.1) then

!********************************** GAUCHE +
    call mrrrr(m,tp,mm)
    if(i.gt.3)Uconf=Uconf+Ug(segpoorsv(i-1,NC))

  elseif (state.eq.2) then
!********************************** GAUCHE -
    call mrrrr(m,tm,mm)
    if(i.gt.3)Uconf=Uconf+Ug(segpoorsv(i-1,NC))

  endif

do j = 1, nbranches(NC) ! loop over branches
 if(branch_pos(j,NC).eq.i-1) then ! last segment was a branching point
   m_branch(:,:,j) = m(:,:) ! save rotation matrix
   xend_branch(:,j) = xend(:,i-1) ! save last position
   state_branch(j) = state
 endif
enddo

m = mm ! update rotation matrix

x(1)=m(1,1)*lseg
x(2)=m(2,1)*lseg
x(3)=m(3,1)*lseg

xend(1,i)=xend(1,i-1)+x(1)   ! ith postion chain
xend(2,i)=xend(2,i-1)+x(2)
xend(3,i)=xend(3,i-1)+x(3)

enddo

!!! Add branches

i = long(NC)-long_branches(NC)

do j = 1, nbranches(NC) 

m(:,:) = m_branch(:,:,j)
xendt(:) = xend_branch(:,j)

do k1=1, branch_long(j,NC)       
i = i + 1

if (k1.eq.1) then ! only first segment of branch
  state = state_branch(j)
  do while (state.eq.state_branch(j)) ! choose a state different to that of the backbone
  rn=rands(seed)
  state=int(rn*3)        ! random select the state= {trans,gauch+,gauch-}
  if (state.eq.3) then
    state=2
  endif
  enddo
else
  rn=rands(seed)
  state=int(rn*3)        ! random select the state= {trans,gauch+,gauch-}
  if (state.eq.3) then
    state=2
  endif
endif

  if (state.eq.0) then
    call mrrrr(m,tt,mm)
  elseif (state.eq.1) then
    call mrrrr(m,tp,mm)
  elseif (state.eq.2) then
    call mrrrr(m,tm,mm)
  endif

do jj = j+1, nbranches(NC) ! loop over remaining branches
 if(branch_pos(jj,NC).eq.i-1) then ! last segment was a branching point
   m_branch(:,:,jj) = m(:,:) ! save rotation matrix
   xend_branch(:,jj) = xend(:,i-1) ! save last position
   state_branch(jj) = state
 endif
enddo

m = mm ! update rotation matrix

x(1)=m(1,1)*lseg
x(2)=m(2,1)*lseg
x(3)=m(3,1)*lseg

xend(1,i)=xendt(1)+x(1)   ! ith postion chain
xend(2,i)=xendt(2)+x(2)
xend(3,i)=xendt(3)+x(3)

xendt(:) = xend(:,i)

enddo ! k1
enddo ! j

dista=0.0                       ! check self avoiding constraint (segmentos)

do k1=1,long(NC)
  do k2=k1+1,long(NC)
    dista=(xend(1,k2)-xend(1,k1))**(2.0)
    dista=dista+(xend(2,k2)-xend(2,k1))**(2.0)
    dista=dista+(xend(3,k2)-xend(3,k1))**(2.0)
    dista=sqrt(dista)+tolerancia
    if (dista.lt.lseg) then
       goto 223
    endif
  enddo
enddo

do i=1,long(NC)
  seglength(segpoorsv(i,NC))=seglength(segpoorsv(i,NC))+1
  do j=1,long(NC)

    if (i.ne.j) then
      distance(i,j)=((xend(1,i)-xend(1,j))**2.0+(xend(2,i)-xend(2,j))**2.0+(xend(3,i)-xend(3,j))**2.0)
      distance(i,j)=sqrt(distance(i,j))
      if (segpoorsv(i,NC).eq.segpoorsv(j,NC))Rgyr(segpoorsv(i,NC))=Rgyr(segpoorsv(i,NC))+distance(i,j)**2.0
      Rgyr(Npoorsv+1)=Rgyr(Npoorsv+1)+distance(i,j)**2.0
      Ugyr=Ugyr-0.5*st(segpoorsv(i,NC),segpoorsv(j,NC))*(lseg/distance(i,j))**(dimf(segpoorsv(i,NC),segpoorsv(j,NC)))
    endif

  enddo
enddo

do is=0,Npoorsv
  Rgyr(is)=sqrt(Rgyr(is)/2.0)
  Rgyr(is)=Rgyr(is)/float(seglength(is))
enddo 

Rgyr(Npoorsv+1)=sqrt(Rgyr(Npoorsv+1)/2.0)
Rgyr(Npoorsv+1)=Rgyr(Npoorsv+1)/float(long(NC))

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
if (ncha.eq.0) goto 223

return
end


subroutine com(xend,xendcom,long)
use seed1
use layer
implicit none

integer i, k, long
real*8 xend(3,long+5),xendcom(3,long+5)
real*8 cm(3)

cm = 0.0

do k = 1,long
  do i = 1,3
    cm(i) = cm(i) + xend(i,k)
  enddo
enddo

cm = cm/float(long)

do k = 1, long
  do i = 1,3
     xendcom(i,k) = xend(i,k) - cm(i)
  enddo
enddo

end subroutine

subroutine print_ent2(xend, indexncha, NC)

use longs
implicit none

real*8 xend(3,200)
integer i,j,jj, indexncha, NC
character*25 filename


! Imprime cadenas en formato ENT

write(filename,'(A6,A1, I3.3, A1, I3.3, A4)') 'cadena','.', indexncha,'.',NC, '.ent'

open(unit=4400, file=filename)

do i=1, long(NC) ! Imprime todo
WRITE(4400,'(A6,I5,A3,I12,A4,F8.3,F8.3,F8.3)') &
"HETATM",i,"  C",i,"    ",xend(1, i)*10,  &
xend(2, i)*10,xend(3, i)*10
end do

i = 1
do j = 1, long(NC)-long_branches(NC)-1 ! Une segmentos backbone
   WRITE((4400),'(A6,I5,I5)')"CONECT", i, i+1
   i = i + 1
end do

do j = 1, nbranches(NC) ! loop over branches
   WRITE((4400),'(A6,I5,I5)')"CONECT", branch_pos(j,NC), i+1
   i = i + 1
   do jj = 1, branch_long(j,NC)-1
       WRITE((4400),'(A6,I5,I5)')"CONECT", i, i+1
       i = i + 1
   enddo
enddo 

WRITE(4400,*)"END"

close(4400)
end subroutine


