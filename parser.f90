!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  parser.f90 = this file reads the inputs from DEFINITIONS.txt
!  
!  FOLLOW THE FORMAT WHEN ADDING NEW INPUTS:
!
!  When implementing new features, please check against test
!  cases and be sure there is backwards compatibility
!
!  Add input validation when possible
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
integer temp


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

! Allocate variables that depend on Ncomp

allocate(flagMD(Ncomp),flagreflex(Ncomp))
allocate(nrot(Ncomp),nrot_corr(Ncomp))
allocate(totalcuantas(Ncomp),cuantas(Ncomp))
allocate(minntotR(Ncomp),minntotZ(Ncomp),maxntotR(Ncomp),maxntotZ(Ncomp))
allocate(long(NComp))
allocate(nbranches(NComp))
allocate(long_branches(NComp))
allocate(npolratio(NComp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! values for not defined variables, change if any variable can take the value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ndi = -10000
ndr = -1.0d10

!!!!!!!!!!! Here all variables are assinged ndi or ndr depending on type

flagMD(:)=ndi
nrot(:)=ndi
flagreflex(:)=ndi
PBCflag = ndi ! flag for PBC in z direction

r_pos = ndr
r_neg = ndr

Npoorsv = ndi ! zero by default
Nacids = ndi
Nbasics = ndi

lseg = ndr
lsegkai = ndr
curvature = ndi

dielP = ndr

dimR = ndi
dimZ = ndi
maxntotR = ndi
maxntotZ = ndi

minntotR=ndi
minntotZ=ndi
dimRini = ndi

deltaR = ndr
deltaZ = ndr

totalcuantas = ndi
long = ndi
npolratio = ndr

flagonekais = ndi
infile = ndi

vtkflag = ndi
maxT = ndi
entflag = ndi

flagkai = ndi

flagtorsionstate = ndi ! zero by defaul

npolini = ndr
npolfirst = ndr
npollast = ndr
npolstep = ndr

Xulimit = ndi
MCfactor = ndi

pHbulk = ndr
csalt = ndr

nbranches = ndi

ta = ndr ! 112


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start reading DEFINITIONS.txt here
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


select case (label)

! flagMD : decides if conformations are read from MD file, see tutorial 
  case('flagMD')
   do NC = 1, NComp
     read(fh, *) flagMD(NC)
   enddo

! nrot : number of rotations
  case('nrot')
   do NC =1,NComp
     read(fh, *) nrot(NC)
   enddo

! flagreflex : reflex conformations (1: no reflex, 2: reflex) 
  case('flagreflex')
   do NC =1,NComp
     read(fh, *) flagreflex(NC)
   enddo

! PBC in z direction, 1: PBC, 2: REFLECTION
  case('PBCflag')
   read(buffer, *, iostat=ios) PBCflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! r_pos, r_neg: salt radius
  case('rsalt')
   read(buffer, *, iostat=ios) r_pos, r_neg
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! Npoorsv : number of Npoorsv types
   case ('Npoorsv')
   read(buffer, *, iostat=ios) Npoorsv
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

   allocate(vpol(0:Npoorsv))
   vpol = ndr

   allocate(Ut(0:Npoorsv),Ug(0:Npoorsv))
   Ut=ndr
   Ug=ndr

   allocate(dimfkais(0:Npoorsv,0:Npoorsv),dimf(0:Npoorsv,0:Npoorsv))
   dimf(:,:) = ndr


! dimf : exponent of vdW interactions
  case ('dimf')
  if(Npoorsv.eq.ndi)call stopparser('Define Npoorsv before reading dimf')
  read(buffer, *, iostat=ios) temp
  if(Npoorsv.ne.temp)call stopparser('Check number of dimf data')
  if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i = 1, Npoorsv
   read(fh,*)(dimf(i,j), j = 1, i)
     do j = 1, i
      dimf(j,i) = dimf(i,j)
     enddo
   enddo

! Utg: trans - gauche energy
   case ('Utg')
   if(Npoorsv.eq.ndi)call stopparser('Define Npoorsv before reading Utg')
   read(buffer, *, iostat=ios) temp
   if(Npoorsv.ne.temp)call stopparser('Check number of Utg data')
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i = 0, Npoorsv
     read(fh,*) Ut(i), Ug(i)
   enddo

! vpol : segment volumens defined according poor-sv
  case ('vpol')
  if(Npoorsv.eq.ndi)call stopparser('Define Npoorsv before reading vpol')
   read(buffer, *, iostat=ios) temp
  if(Npoorsv.ne.temp)call stopparser('Check number of vpol data')
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i=0,npoorsv
     read(fh, *) vpol(i)
   enddo


! Nacids: number of acid molecules and pKas
  case ('Nacids')
   read(buffer, *, iostat=ios) Nacids
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(pKa(Nacids),Ka(Nacids))

   allocate(vpol_a(Nacids))
   vpol_a(:) = ndr

   do i=1,Nacids
     read(fh,*)pKa(i) ! acid constants of each acid segment
   enddo

! Nbasics: number of basic molecules and pKas
  case ('Nbasics')
   read(buffer, *, iostat=ios) Nbasics
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   allocate(Kb(Nbasics),pKb(Nbasics))


   allocate(vpol_b(Nbasics))
   vpol_b(:) = ndr

   do i=1,Nbasics
     read(fh,*)pKb(i) ! acid constants of each acid segment
   enddo

 
! vpol_a : segment volumens defined according acid-base
  case('vpol_a')
  if(Nacids.eq.ndi)call stopparser('Define Nacids before reading vpol_a')
   read(buffer, *, iostat=ios) temp
  if(Nacids.ne.temp)call stopparser('Check number of vpol_a data')
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i=1,Nacids
     read(fh,*) vpol_a(i)
   enddo

! vpol_b : segment volumens defined according acid-base
  case('vpol_b')
  if(Nbasics.eq.ndi)call stopparser('Define Nbasics before reading vpol_a')
   read(buffer, *, iostat=ios) temp
  if(Nbasics.ne.temp)call stopparser('Check number of vpol_b data')
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   do i=1,Nbasics
     read(fh,*) vpol_b(i)
   enddo

! lseg : segment length for chain generation
  case ('lseg')
   read(buffer, *, iostat=ios) lseg
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! lsegkai : segment length for chi calculations
  case ('lsegkai')
   read(buffer, *, iostat=ios) lsegkai
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)


! curvature
! 0 = plane, 1: cylinder, 2 : micelle, 3 : plane with PBC
  case ('curvature')
   read(buffer, *, iostat=ios) curvature
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! Dielectric constant of the polymer
  case ('dielP')
   read(buffer, *, iostat=ios) dielP
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! dimensions, sizes in R, Z and standard maxntot (maximum values of R and Z for the CM)
  case ('dimensions')
   read(buffer, *, iostat=ios) dimR, dimZ, maxntotR_all, maxntotZ_all
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)
   ntot=dimR*dimZ
   maxntotR(:)=maxntotR_all
   maxntotZ(:)=maxntotZ_all
   maxntotR_max=maxntotR_all

! mintot (minimum values of R and Z for the CM)
  case('minntot')
  if((dimR.eq.ndi).or.(dimZ.eq.ndi))call stopparser('Define dimR and dimZ before minntot')
    do NC=1,Ncomp
      read(fh, *) minntotR(NC),minntotZ(NC)
    enddo
    minntotR_min = minval(minntotR)

 ! maxntot (maximum values of R and Z for the CM)
  case('maxntot') ! must be read after reading the lattice properties (case dimensions)
  if((dimR.eq.ndi).or.(dimZ.eq.ndi))call stopparser('Define dimR and dimZ before maxntot')
    do NC=1,Ncomp
      read(fh, *) maxntotR(NC),minntotZ(NC)
    enddo
    maxntotR_max = maxval(maxntotR)

! dimRini : initial value of R, useful for curved lamella.   
  case ('dimRini')
   read(buffer, *, iostat=ios) dimRini
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! layersize, deltaR and deltaZ : size of the layers in nm units
  case ('layersize')
   read(buffer, *, iostat=ios) deltaR, deltaZ
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),'= ',trim(buffer)

! cuantas: number of conformations
  case ('cuantas')
    do NC = 1, NComp
      read(fh, *) totalcuantas(NC)
    enddo

! long : length of chains
  case ('long')
   do NC = 1, NComp
     read(fh, *) long(NC)
   enddo
   maxlong = maxval(long)

! npolratio : multiplicative factor for the number of molecules of type NC  
  case ('npolratio')
   do NC = 1, NComp
     read(fh, *) npolratio(NC)
   enddo

! flagonekais: read all kais from kai.001.001.dat   
 case ('flagonekais')
   read(buffer, *, iostat=ios) flagonekais
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! infile : 0: do not use initial guess, 1: use initial guess from fort files, 2: use initial guess from in.in
 case ('infile')
   read(buffer, *, iostat=ios) infile
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! vtkflag : save results as vtk?   
  case ('vtkflag')
   read(buffer, *, iostat=ios) vtkflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! entflag : save chain conformations as *.ent files for visualization
  case ('entflag')
   read(buffer, *, iostat=ios) entflag
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! maxT : number of rotations in vtk file, useful to make nice figures   
  case ('vtkT')
   read(buffer, *, iostat=ios) maxT
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! flagkai: 0: read kai from files, 1: generate new kai tables
  case ('flagkai')
   read(buffer, *, iostat=ios) flagkai

! flagtorsionstate: 0: random dihedral, 1: read dihedrals from file 
  case ('flagtorsionstate')
   read(buffer, *, iostat=ios) flagtorsionstate

! npol: Aggregtion number or density
  case ('npol')
   read(buffer, *, iostat=ios) npolini, npolfirst, npollast, npolstep
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! Xulimit : cutoff of vdW integration   
  case ('Xulimit')
   read(buffer, *, iostat=ios) Xulimit
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! MCsteps : number of steps in vdW integration   
  case('MCsteps')
   read(buffer, *, iostat=ios) MCfactor
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! csalt : salt concentration   
  case ('csalt')
   read(buffer, *, iostat=ios) csalt
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! pHbulk : bulk pH   
  case ('pHbulk')
   read(buffer, *, iostat=ios) pHbulk
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

! number and position of branches   
   case ('nbranches')
  
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

! torsion angle   
  case('torsion_angle')
   read(buffer, *, iostat=ios) ta
   if(rank.eq.0)write(stdout,*) 'parser:','Set ',trim(label),' = ',trim(buffer)

endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check validity of input
!  
! If critical variables are undefined, then stop
! 
! Otherwise, assign default values
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! flagMD
do NC=1,Ncomp
  if(flagMD(NC).eq.ndi) then
     flagMD(NC) = 0
     if(rank.eq.0)write(stdout,*) 'Component ', NC, ' flagMD undefined, use default value (0 : do not read MD conformations)'
   endif
enddo          
        
! nrot
do NC=1,Ncomp
  if(nrot(NC).eq.ndi) then
     nrot(NC) = 12
     if(rank.eq.0)write(stdout,*) 'Component ', NC, ' nrot undefined, use default value (12)'
   endif
enddo          

! flagreflex
do NC=1,Ncomp
  if(flagreflex(NC).eq.ndi) then
     flagreflex(NC) = 1
     if(rank.eq.0)write(stdout,*) 'Component ', NC, ' flagreflex undefined, use default value (1 : do not reflex conformations)'
   endif
enddo          

! PBCflag
if(PBCflag.eq.ndi) then
  PBCflag = 1        
  if(rank.eq.0)write(stdout,*) 'PBCflag undefined, use default value (1 : PBC in z)'
endif

! r_pos        
if(r_pos.eq.ndr) then
  r_pos = 0.3  
  if(rank.eq.0)write(stdout,*) 'r_pos undefined, use default value (0.3)'
endif

! r_neg       
if(r_neg.eq.ndr) then
  r_neg = 0.3  
  if(rank.eq.0)write(stdout,*) 'r_neg undefined, use default value (0.3)'
endif

! Npoorsv
if(Npoorsv.eq.ndi)call stopundef('Npoorsv')

! Nacids
if(Nacids.eq.ndi)call stopundef('Nacids')

! Nbasics
if(Nbasics.eq.ndi)call stopundef('Nbasics')

! dimf
do i=0,Npoorsv
  do j=0,Npoorsv
     if(dimf(i,j).eq.ndr) then
        dimf(i,j) = 6.0
        if(rank.eq.0)write(stdout,*) 'dimf ', i,j, ' undefined, use default value (6.0)'
     endif
  enddo          
enddo          

! Ut / Ug
do i=0,Npoorsv
     if(Ut(i).eq.ndr) then
        Ut(i) = 0.0
        if(rank.eq.0)write(stdout,*) 'Ut ', i,' undefined, use default value (0.0)'
     endif
     if(Ug(i).eq.ndr) then
        Ug(i) = 0.0
        if(rank.eq.0)write(stdout,*) 'Ug ', i,' undefined, use default value (0.0)'
     endif
enddo          

! vpol
do i=0,Npoorsv
  if(vpol(i).eq.ndr)call stopundef('vpol')
enddo          

! vpol_a
do i=1,Nacids
  if(vpol_a(i).eq.ndr)call stopundef('vpol_a')
enddo          

! vpol_b
do i=1,Nbasics
  if(vpol_b(i).eq.ndr)call stopundef('vpol_b')
enddo          

! lseg
if(lseg.eq.ndr)call stopundef('lseg')

! lsegkai
if(lsegkai.eq.ndr)call stopundef('lsegkai')

! curvature
if(curvature.eq.ndi)call stopundef('curvature')
select case (curvature)
   case(0,1,2)
   case DEFAULT
       call stopparser('Curvature should be 0,1,2 or 3') 
endselect

! dielP   
if(dielP.eq.ndr) then
  dielP = 78.54
  if(rank.eq.0)write(stdout,*) 'dielP undefined, use default value (78.54)'
endif

! dimR
if(dimR.eq.ndi)call stopundef('dimR')

! dimZ
if(dimZ.eq.ndi)call stopundef('dimZ')

! maxntotR, maxntotZ
do NC=1,Ncomp
  if(maxntotR(NC).eq.ndi)call stopundef('maxntotR')
  if(maxntotZ(NC).eq.ndi)call stopundef('maxntotZ')
  if(totalcuantas(NC).eq.ndi)call stopundef('cuantas')
enddo

! minntotR
do NC=1,Ncomp
  if(minntotR(NC).eq.ndi) then
     minntotR(NC) = 1
     if(rank.eq.0)write(stdout,*) 'Component ', NC, ' minntotR undefined, use default value (1)'
   endif
enddo          
minntotR_min = minval(minntotR)

! minntotZ
do NC=1,Ncomp
  if(minntotZ(NC).eq.ndi) then
     minntotZ(NC) = 1
     if(rank.eq.0)write(stdout,*) 'Component ', NC, ' minntotZ undefined, use default value (1)'
   endif
enddo          

! dimRini
if(dimRini.eq.ndi) then
  dimRini = 0
  if(rank.eq.0)write(stdout,*) 'dimRini undefined, use default value (0)'
endif

! deltaR and deltaZ
if(deltaR.eq.ndr)call stopundef('deltaR')
if(deltaZ.eq.ndr)call stopundef('deltaZ')

! totalcuantas
do NC=1,Ncomp
  if(totalcuantas(NC).eq.ndi)call stopundef('cuantas')
enddo          

! long
do NC=1,Ncomp
  if(long(NC).eq.ndi)call stopundef('long')
enddo          

! npolration
do NC=1,Ncomp
  if(npolratio(NC).eq.ndr)call stopundef('npolratio')
enddo          

! flagonekais
if(flagonekais.eq.ndi) then
  flagonekais = 0
  if(rank.eq.0)write(stdout,*) 'flagonekais undefined, use default value (0 : read kais from different files kai.*.*.dat)'
endif

! infile
if(infile.eq.ndi)call stopundef('infile')
select case (infile)
   case(0,1,2)
   case DEFAULT
       call stopparser('infile should be 0,1 or 2') 
endselect

! vtkflag
if(vtkflag.eq.ndi) then
  vtkflag = 0
  if(rank.eq.0)write(stdout,*) 'vtkflag undefined, use default value (0 : do not save vtk)'
endif

! maxT
if(maxT.eq.ndi) then
  maxT = 1
  if(rank.eq.0)write(stdout,*) 'vtkT undefined, use default value (1)'
endif

! entflag
if(entflag.eq.ndi) then
  entflag = 0
  if(rank.eq.0)write(stdout,*) 'entflag undefined, use default value (0 : do not save conformations as ent files)'
endif

! flagkai
if(flagkai.eq.ndi) then
  flagkai = 0
  if(rank.eq.0)write(stdout,*) 'flagkai undefined, use default value (0 : read kais from file, do not generate)'
endif

! flagtorsionstate
if(flagtorsionstate.eq.ndi) then
  flagtorsionstate = 0
  if(rank.eq.0)write(stdout,*) 'flagtorsionstate undefined, use default value (0 : do not read torsion states, use random)'
endif

! npol
if(npolini.eq.ndr)call stopundef('npolini')
if(npolfirst.eq.ndr)call stopundef('npolfirst')
if(npollast.eq.ndr)call stopundef('npollast')
if(npolstep.eq.ndr)call stopundef('npolstep')

! Xulimit
if(Xulimit.eq.ndi)call stopundef('Xulimit')

! MCfactor/MCsteps
if(MCfactor.eq.ndi) then
  MCfactor = 60
  if(rank.eq.0)write(stdout,*) 'MCsteps undefined, use default value (60)'
endif


! pH
if(pHbulk.eq.ndr)call stopundef('pHbulk')

! Csalt
if(csalt.eq.ndr)call stopundef('csalt')

! nbranches
do NC=1,Ncomp
  if(nbranches(NC).eq.ndi) then
     nbranches(NC) = 0
     if(rank.eq.0)write(stdout,*) 'Component ', NC, ' nbranches undefined, use default value (0)'
   endif
enddo          

! ta
if(ta.eq.ndr) then
  ta = 112.0
  if(rank.eq.0)write(stdout,*) 'torsion-angle undefined, use default value (112.0)'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                                        Auxiliary calculations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Correct number of rotations if reflection is used

do NC=1,Ncomp
  nrot_corr(NC)=nrot(NC)*flagreflex(NC)
enddo


! block data for MPI implementation
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                             Read chain structure from structure.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     Read st(i,j) from epsilon.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

subroutine stopparser(namevar)
character(len=*) :: namevar
integer, parameter :: stdout = 6
write(stdout,*) 'parser:', namevar
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
end

