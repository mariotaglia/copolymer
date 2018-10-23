subroutine parser
use mcharge
use longs
use globals
use bulk
use MPI
use mkai
use volume

implicit none
integer block_cuantas, restcuantas

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer, parameter :: stdout = 6
integer ios
integer line, linemax
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
character basura
integer ndi ! undetermined integer
real*8 ndr ! undetermined real

! not defined variables, change if any variable can take the value
ndi = -1e5
ndr = -1.0d10

! default values, if ndi or ndr is used, then variable is required
dielP = 78.54
curvature = ndi
ntot = ndi
maxntotcounter_ini = ndi
maxntot = ndi
totalcuantas = ndi
long = ndi
Npoorsv = ndi ! zero by default
Nacids = 0
Nbasics = 0
infile = ndi
flagkai = 0 ! zero by default
npolini = ndi
npolfirst = ndi
npollast = ndi
npolstep = ndi
Xulimit = ndi
lseg = ndr
Csalt = ndr
pHbulk = ndr
nbranches = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)write(stdout,*) 'parser:', 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


select case (label)

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

  case ('dielP')
   read(buffer, *, iostat=ios) dielP
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('curvature')
   read(buffer, *, iostat=ios) curvature
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('ntot')
   read(buffer, *, iostat=ios) ntot, maxntotcounter_ini, maxntot
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('cuantas')
   read(buffer, *, iostat=ios) totalcuantas
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
  ! Auxiliary calculations
   block_cuantas=int(totalcuantas/size/12)
   cuantas=block_cuantas*12
   restcuantas=totalcuantas-size*12*block_cuantas
   if(rank.eq.(size-1))cuantas=cuantas+restcuantas

  case ('long')
   read(buffer, *, iostat=ios) long
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

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

  case ('infile')
   read(buffer, *, iostat=ios) infile
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
     read(fh,*), Ut(i), Ug(i)
   enddo

  case ('csalt')
   read(buffer, *, iostat=ios) Csalt
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('pHbulk')
   read(buffer, *, iostat=ios) pHbulk
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

  case ('nbranches')
   read(buffer, *, iostat=ios) nbranches
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   allocate (branch_pos(nbranches))
   allocate (branch_long(nbranches))
   long_branches = 0

   do j = 1, nbranches
   read(fh,*) branch_pos(j), branch_long(j)
   long_branches = long_branches + branch_long(j)
   enddo
endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check validity of input
! 

if(curvature.eq.ndr)call stopundef('curvature')
if(ntot.eq.ndi)call stopundef('ntot')
if(maxntotcounter.eq.ndi)call stopundef('maxntot_counter')
if(maxntot.eq.ndi)call stopundef('maxntot')
if(totalcuantas.eq.ndi)call stopundef('cuantas')
if(long.eq.ndi)call stopundef('long')
if(infile.eq.ndi)call stopundef('infile')
if(npolini.eq.ndi)call stopundef('npolini')
if(npolfirst.eq.ndi)call stopundef('npolfirst')
if(npollast.eq.ndi)call stopundef('npollast')
if(npolstep.eq.ndi)call stopundef('npolstep')
if(Xulimit.eq.ndi)call stopundef('Xulimit')
if(lseg.eq.ndr)call stopundef('lseg')
if(pHbulk.eq.ndi)call stopundef('pHbulk')

if (rank.eq.0) then
 print*, 'pKa = ', pka
 print*, 'pKb = ', pKb
 print*, 'vpol = ', vpol
 print*, 'vpol_a = ', vpol_a
 print*, 'vpol_b = ', vpol_b
 print*, 'dimf = ', dimf
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read chain structure from structure.in

open(file='structure.in', unit = 9)
allocate(segpoorsv(long))
allocate(acidtype(long))
allocate(basictype(long))
do i = 1, long
read(9,*),segpoorsv(i), acidtype(i), basictype(i)
enddo

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

