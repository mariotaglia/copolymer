subroutine parser
use mcharge
use longs
use globals
use bulk
use MPI
use mkai
use volume
use layer
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


print*, "LEO ARCHVIO" !LEO
! not defined variables, change if any variable can take the value
ndi = -10000
ndr = -1.0d10

! default values, if ndi or ndr is used, then variable is required
dielP = 78.54
curvature = ndi
ntot = ndi
maxntotR = ndi
maxntotZ = ndi
dimR = ndi
dimZ = ndi
totalcuantas = ndi
Npoorsv = ndi ! zero by default
Nacids = 0
Nbasics = 0
pKcopmol = ndr ! LEO
pKcopion = ndr ! LEO
pKmolion = ndr ! LEO
infile = ndi
flagkai = 0 ! zero by default
r_pos = 0.3
r_neg = 0.3
npolini = ndi
npolfirst = ndi
npollast = ndi
npolstep = ndi
Xulimit = ndi
lseg = ndr
Csalt = ndr
pHbulk = ndr
PBCflag = 1 ! flag for PBC in z direction
vtkflag = 0
entflag = 0
maxT = 1

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

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


select case (label)

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
  
  case ('vcopmol')
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   read(fh, *) vcopmol

  case ('lseg')
   read(buffer, *, iostat=ios) lseg
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
  
  case ('layersize')
   read(buffer, *, iostat=ios) deltaR, deltaZ
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),'= ',trim(buffer)

  case ('cuantas')
   read(buffer, *, iostat=ios) totalcuantas
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  ! Auxiliary calculations
   block_cuantas=int(totalcuantas/size/12)
   cuantas=block_cuantas*12
   restcuantas=totalcuantas-size*12*block_cuantas
   if(rank.eq.(size-1))cuantas=cuantas+restcuantas

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
  
  case ('npoorsv')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  
  case ('dimf')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(dimfkais(0:Npoorsv,0:Npoorsv),dimf(0:Npoorsv,0:Npoorsv))
   dimf(0,0)=0
   do i = 1, Npoorsv
   dimf(0,i)=0.
   dimf(i,0)=0.
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

  case ('pKcopmol') ! LEO
   read(buffer, *, iostat=ios) pKcopmol
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('pKcopion') ! LEO
   read(buffer, *, iostat=ios) pKcopion
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  
  case ('pKmolion') ! LEO
   read(buffer, *, iostat=ios) pKmolion
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   
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

  case ('npol')
   read(buffer, *, iostat=ios) npolini, npolfirst, npollast, npolstep
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('Xulimit')
   read(buffer, *, iostat=ios) Xulimit
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('Utg')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(Ut(0:npoorsv),Ug(0:Npoorsv))
   Ut=0.0
   Ug=0.0
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
if(totalcuantas.eq.ndi)call stopundef('cuantas')
if(infile.eq.ndi)call stopundef('infile')
if(npolini.eq.ndi)call stopundef('npolini')
if(npolfirst.eq.ndi)call stopundef('npolfirst')
if(npollast.eq.ndi)call stopundef('npollast')
if(npolstep.eq.ndi)call stopundef('npolstep')
if(Xulimit.eq.ndi)call stopundef('Xulimit')
if(lseg.eq.ndr)call stopundef('lseg')
if(pHbulk.eq.ndi)call stopundef('pHbulk')
if(pKcopmol.eq.ndr)call stopundef('pKcopmol') ! LEO
if(pKcopion.eq.ndr)call stopundef('pKcopion') ! LEO
if(pKmolion.eq.ndr)call stopundef('pKmolion') ! LEO

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

do NC = 1, Ncomp

write(filename2,'(A10,I3.3,A3)')'structure.',NC,'.in'
open(file=filename2, unit = 9)
do i = 1, long(NC)
read(9,*)segpoorsv(i,NC), acidtype(i,NC), basictype(i,NC)
enddo
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

