subroutine parser
use mcharge
use longs
use globals
use bulk
use MPI
use mkai
use volume
use layer
use cadenaMD
use senos
implicit none
integer block_cuantas, restcuantas

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer, parameter :: stdout = 6
integer ios
integer line
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
integer ndi ! undetermined integer
real*8 ndr ! undetermined real
integer NC
character*16 filename2

! not defined variables, change if any variable can take the value
ndi = -10000
ndr = -1.0d10

! default values, if ndi or ndr is used, then variable is required
dielP = 78.54
curvature = ndi
ntot = ndi
minntotR = 1
minntotZ = 1
maxntotR = ndi
maxntotZ = ndi
dimR = ndi
dimRini = 0
dimZ = ndi
!totalcuantas = ndi
Npoorsv = ndi ! zero by default
Nacids = 0
Nbasics = 0
infile = ndi
flagkai = 0 ! zero by default
flagtorsionstate = 0 ! zero by defaul
r_pos = 0.3
r_neg = 0.3
MCfactor = 60
npolini = ndi
npolfirst = ndi
npollast = ndi
npolstep = ndi
Xulimit = ndi
lseg = ndr
lsegkai = ndr
Csalt = ndr
pHbulk = ndr
PBCflag = 1 ! flag for PBC in z direction
vtkflag = 0
entflag = 0
maxT = 1
!flagMD = 0
!nrot = 12 ! default number of rotations
ta = 112

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)write(stdout,*) 'parser:', 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

read(fh, '(A)', iostat=ios) buffer
pos = scan(buffer, ' ')
label = buffer(1:pos)
buffer = buffer(pos+1:)

if(label.ne.'Ncomp') then
   if(rank.eq.0)print*,'FIRST LINE IN DEFINITIONS SHOULD BE Ncomp'
   stop
endif

read(buffer, *, iostat=ios) Ncomp


allocate(flagMD(Ncomp),flagreflex(Ncomp))
allocate(nrot(Ncomp),nrot_corr(Ncomp))
allocate(totalcuantas(Ncomp),cuantas(Ncomp))
nrot(:)=12
flagMD(:)=0
flagreflex(:)=1
totalcuantas(:)=ndi




do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


select case (label)
 
  case('flagMD')
   do NC = 1, NComp
     read(fh, *) flagMD(NC)
   enddo
  
  case('nrot')
   do NC =1,NComp
     read(fh, *) nrot(NC)
   enddo
 
  case('flagreflex')
   do NC =1,NComp
     read(fh, *) flagreflex(NC)
   enddo

  case('PBCflag')
   read(buffer, *, iostat=ios) PBCflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case('rsalt')
   read(buffer, *, iostat=ios) r_pos, r_neg
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('vpol')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(vpol(0:Npoorsv))
   do i=0,npoorsv
     read(fh, *) vpol(i)
   enddo

  case('vpol_a')
   read(buffer, *, iostat=ios) Nacids
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(vpol_a(Nacids))
   do i=1,Nacids
     read(fh,*) vpol_a(i)
   enddo

  case('vpol_b')
   read(buffer, *, iostat=ios) Nbasics
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(vpol_b(Nbasics))
   do i=1,Nbasics
     read(fh,*) vpol_b(i)
   enddo

  case ('lseg')
   read(buffer, *, iostat=ios) lseg
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('lsegkai')
   read(buffer, *, iostat=ios) lsegkai
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('dielP')
   read(buffer, *, iostat=ios) dielP
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('curvature')
   read(buffer, *, iostat=ios) curvature
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('dimensions')
   read(buffer, *, iostat=ios) dimR, dimZ, maxntotR, maxntotZ
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   ntot=dimR*dimZ

  case('minntot')
   read(buffer, *, iostat=ios) minntotR, minntotZ
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('dimRini')
   read(buffer, *, iostat=ios) dimRini
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


  case ('layersize')
   read(buffer, *, iostat=ios) deltaR, deltaZ
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),'= ',trim(buffer)

  case ('cuantas')
    do NC = 1, NComp
      read(fh, *) totalcuantas(NC)
    enddo

  case ('long')
   allocate(long(NComp))
   do NC = 1, NComp
     read(fh, *) long(NC)
   enddo
   maxlong = maxval(long)

  case ('npolratio')
   allocate(npolratio(NComp))
   do NC = 1, NComp
     read(fh, *) npolratio(NC)
   enddo
  
  case ('Npoorsv')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(Ut(0:Npoorsv),Ug(0:Npoorsv))
   Ut=0.0
   Ug=0.0
   allocate(dimfkais(0:Npoorsv,0:Npoorsv),dimf(0:Npoorsv,0:Npoorsv))
   dimf(:,:)=6.
   dimf(0,0)=0.
   print*,dimf

  case ('dimf')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i = 1, Npoorsv
   read(fh,*)(dimf(i,j), j = 1, i)
     do j = 1, i
      dimf(j,i) = dimf(i,j)
     enddo
   enddo

  case ('Nacids')
   read(buffer, *, iostat=ios) Nacids
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(pKa(Nacids),Ka(Nacids))
   do i=1,Nacids
     read(fh,*)pKa(i) ! acid constants of each acid segment
   enddo

  case ('Nbasics')
   read(buffer, *, iostat=ios) Nbasics
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(Kb(Nbasics),pKb(Nbasics))
   do i=1,Nbasics
     read(fh,*)pKb(i) ! acid constants of each acid segment
   enddo

  case ('infile')
   read(buffer, *, iostat=ios) infile
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('vtkflag')
   read(buffer, *, iostat=ios) vtkflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('entflag')
   read(buffer, *, iostat=ios) entflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('vtkT')
   read(buffer, *, iostat=ios) maxT
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


  case ('flagkai')
   read(buffer, *, iostat=ios) flagkai

  case ('flagtorsionstate')
   read(buffer, *, iostat=ios) flagtorsionstate

  case ('npol')
   read(buffer, *, iostat=ios) npolini, npolfirst, npollast, npolstep
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('Xulimit')
   read(buffer, *, iostat=ios) Xulimit
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  
  case('MCsteps')
   read(buffer, *, iostat=ios) MCfactor
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('Utg')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i = 0, Npoorsv
     read(fh,*) Ut(i), Ug(i)
   enddo

  case ('csalt')
   read(buffer, *, iostat=ios) Csalt
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('pHbulk')
   read(buffer, *, iostat=ios) pHbulk
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   case ('nbranches')
   allocate(nbranches(NComp))
   allocate(long_branches(NComp))
   
   do NC = 1, NComp
     read(fh, *)nbranches(NC)
   enddo
   maxnbranches = maxval(nbranches)
 
   allocate (branch_pos(maxnbranches,Ncomp))
   allocate (branch_long(maxnbranches,Ncomp))

   do NC = 1, Ncomp
   long_branches(NC) = 0
     do j = 1, nbranches(NC)
        read(fh,*) branch_pos(j,NC), branch_long(j,NC)
        long_branches(NC) = long_branches(NC) + branch_long(j,NC)
     enddo
   enddo ! NC
   
  case('torsion_angle')
   read(buffer, *, iostat=ios) ta
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check validity of input
! 

if(curvature.eq.ndr)call stopundef('curvature')
if(ntot.eq.ndi)call stopundef('ntot')
if(maxntotR.eq.ndi)call stopundef('maxntotR')
if(maxntotZ.eq.ndi)call stopundef('maxntotZ')

do NC=1,Ncomp
  if(totalcuantas(NC).eq.ndi)call stopundef('cuantas')
enddo

if(infile.eq.ndi)call stopundef('infile')
if(npolini.eq.ndi)call stopundef('npolini')
if(npolfirst.eq.ndi)call stopundef('npolfirst')
if(npollast.eq.ndi)call stopundef('npollast')
if(npolstep.eq.ndi)call stopundef('npolstep')
if(Xulimit.eq.ndi)call stopundef('Xulimit')
if(lseg.eq.ndr)call stopundef('lseg')
if(lsegkai.eq.ndr)call stopundef('lsegkai')
if(pHbulk.eq.ndi)call stopundef('pHbulk')


! Auxiliary calculations

do NC=1,Ncomp
  nrot_corr(NC)=nrot(NC)*flagreflex(NC)
enddo
ncha_max=nrot_corr(1)
block_cuantas=int(totalcuantas(1)/size/nrot_corr(1))
cuantas(1)=block_cuantas*nrot_corr(1)
restcuantas=totalcuantas(1)-size*nrot_corr(1)*block_cuantas
if(rank.eq.(size-1))cuantas(1)=cuantas(1)+restcuantas

cuantas_max=cuantas(1)

  do NC=2,Ncomp
     block_cuantas=int(totalcuantas(NC)/size/nrot_corr(NC))
     cuantas(NC)=block_cuantas*nrot_corr(NC)
     restcuantas=totalcuantas(NC)-size*nrot_corr(NC)*block_cuantas
     if(rank.eq.(size-1))cuantas(NC)=cuantas(NC)+restcuantas
     if(cuantas(NC).gt.cuantas(1))cuantas_max=cuantas(NC)
     if(nrot_corr(NC).gt.ncha_max)ncha_max=nrot_corr(NC)
  enddo


! write output

if (rank.eq.0) then
 print*, 'pKa = ', pka
 print*, 'pKb = ', pKb
 print*, 'vpol = ', vpol
 print*, 'vpol_a = ', vpol_a
 print*, 'vpol_b = ', vpol_b
 print*, 'r_pos = ', r_pos
 print*, 'r_neg = ', r_neg
 print*, 'dimf = ', dimf
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read chain structure from structure.in

allocate(segpoorsv(maxlong,Ncomp))
allocate(acidtype(maxlong,Ncomp))
allocate(basictype(maxlong,Ncomp))
allocate(torsionstate(maxlong,Ncomp))

do NC = 1, Ncomp

  write(filename2,'(A10,I3.3,A3)')'structure.',NC,'.in'
  open(file=filename2, unit = 9)

  if (flagMD(NC).eq.0) then ! RIS conformation
    if (flagtorsionstate.eq.0) then
      do i = 1, long(NC)
        read(9,*)segpoorsv(i,NC), acidtype(i,NC), basictype(i,NC) ! , torsionstate(i,NC)
        torsionstate(i,NC)=3
      enddo
    else 
      do i = 1, long(NC)
        read(9,*)segpoorsv(i,NC), acidtype(i,NC), basictype(i,NC), torsionstate(i,NC)
      enddo
    endif

  else
    read(9,*) nMD ! number of MD beads
    do i = 1, nMD
      read(9,*) j, MDHs(i,NC), MDsegpoorsv(i,NC), MDacidtype(i,NC), MDbasictype(i,NC) ! properties of MD bead i, first column is 0 for H or 1 for heavy
    enddo
  endif

  close(9)

enddo ! NC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read st(i,j) from epsilon.in

open(file='epsilon.in', unit=10)
allocate(st(0:Npoorsv,0:Npoorsv))
st(0,0)=0.
do i = 1, Npoorsv
st(0,i)=0.
st(i,0)=0.
read(10,*)(st(i,j), j = 1, i)
 do j = 1, i
   st(j,i) = st(i,j)
 enddo
enddo







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

subroutine stopundef(namevar)
character(len=*) :: namevar
integer, parameter :: stdout = 6
write(stdout,*) 'parser:', 'Variable ', namevar, ' is undefined '
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
end

