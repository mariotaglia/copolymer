!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 1D-coordenadas
! polares usando un metodo de MC
!
!#####################################################################


subroutine kai

use globals
use layer
use volume 
use MPI
use mkai
implicit none
real*8 xmin,xmax,ymin,ymax,zmin,zmax
integer MCsteps ! numero de steps de MC


real*8 R,z
integer ii, is, js, a, b, c
real*8 pi
real*8 x1,x2,y1, y2, z1, z2, vect
integer ix,iy,iz
integer jR, jZ

real*8 cutoff
real*8, allocatable :: sumaXu(:,:)
character*16 kaisfilename
integer, external :: PBCSYMI


if(rank.eq.0)print*,'Kai calculation'


allocate(sumaXu(Npoorsv,Npoorsv))

cutoff = (float(Xulimit)+0.5)*deltaR

pi=dacos(-1.0d0)          ! pi = arccos(-1) 

Xu = 0.0 ! vector Xu

MCsteps = MCfactor*Xulimit
sumaXu(:,:)=0.0

if (flagkai.eq.1) then

select case (curvature)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(0,1,2)
         
   do ii = 1, dimR ! loop sobre cada posicion del segmento

      ymax = cutoff
      ymin = -cutoff
      zmax = cutoff
      zmin = -cutoff

      xmax = cutoff + (dfloat(ii+dimRini) - 0.5)*deltaR
      xmin = -cutoff + (dfloat(ii+dimRini) - 0.5)*deltaR
      
      do ix = 1, MCsteps
      do iy = 1, MCsteps
      do iz = 1, MCsteps

  ! coordenadas del segmento (x1,y1,z1) y del punto a integrar (x2,y2,z2)

         x1 = (dfloat(ii+dimRini) - 0.5)*deltaR ! asume theta segmento = 0, z segmento = 0 y segmento en el centro de la layer
         y1 = 0.0
         z1 = 0.0

         x2 = xmin + (xmax-xmin)*dfloat(ix-1)/dfloat(MCsteps-1)
         y2 = ymin + (ymax-ymin)*dfloat(iy-1)/dfloat(MCsteps-1)
         z2 = zmin + (zmax-zmin)*dfloat(iz-1)/dfloat(MCsteps-1)

         select case (abs(curvature))
         case (0)
         R = abs(x2)
         Z = z2
         case (1)
         R = sqrt(x2**2 + y2**2)
         Z = z2
         case (2)
         R = sqrt(x2**2 + y2**2 + z2**2)
         if(dimZ.ne.1)stop
         Z = 0.0
         end select

         vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diferencia
         jR = int(R/deltaR)-dimRini+1 ! jR tiene la celda donde cae el punto a integrar
         jZ = int(anint(Z/deltaZ)) ! OJO

         if((jR.ge.1).and.(jR.le.dimR)) then
         
         if(vect.le.(cutoff)) then ! esta dentro de la esfera del cut-off   
         

         if(vect.ge.lsegkai) then ! esta dentro de la esfera del segmento

           do is=1,Npoorsv
           do js=1,Npoorsv
              Xu(ii, jR, jZ, is, js) = Xu(ii, jR, jZ, is, js) + ((lsegkai/vect)**dimf(is, js)) ! incluye el jacobiano R(segmento)
           enddo
           enddo
         endif
         endif

         endif
         
      enddo!iz
      enddo!iy
      enddo!ix
      
      do jR = 1, dimR
      do jZ = -Xulimit, Xulimit
         do is=1,Npoorsv
         do js=1,Npoorsv
           Xu(ii, jR, jZ, is, js) = Xu(ii, jR, jZ, is, js)/(MCsteps**3)*(2.0*cutoff)**3
         enddo
         enddo
      enddo
      enddo
  end do ! ii

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(3) ! lamella con PBC
         
   do ii = 1, dimR ! loop over the initial position of each segment

      ymax = cutoff
      ymin = -cutoff
      zmax = cutoff
      zmin = -cutoff

      xmax = cutoff + (dfloat(ii+dimRini) - 0.5)*deltaR
      xmin = -cutoff + (dfloat(ii+dimRini) - 0.5)*deltaR
      
      do ix = 1, MCsteps
      do iy = 1, MCsteps
      do iz = 1, MCsteps

  ! coordenadas del segmento (x1,y1,z1) y del punto a integrar (x2,y2,z2)

  x1 = dfloat(ii)*deltaR ! for ii=1, x1=delta and if x2 is between delta/2 and 3/2*delta, the point falls in layer 1
         y1 = 0.0
         z1 = 0.0

         x2 = xmin + (xmax-xmin)*dfloat(ix-1)/dfloat(MCsteps-1)
         y2 = ymin + (ymax-ymin)*dfloat(iy-1)/dfloat(MCsteps-1)
         z2 = zmin + (zmax-zmin)*dfloat(iz-1)/dfloat(MCsteps-1)

         R = x2
         Z = z2

         vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diferencia

         jZ = int(anint(Z/deltaZ)) 
         jR = int(anint(R/deltaR))

         jR = PBCSYMI(jR,dimR) !  puts the segment within the calculation box using PBC,

         if(vect.le.(cutoff)) then ! esta dentro de la esfera del cut-off   
         if(vect.ge.lsegkai) then ! esta dentro de la esfera del segmento
           do is=1,Npoorsv
           do js=1,Npoorsv
              Xu(ii, jR, jZ, is, js) = Xu(ii, jR, jZ, is, js) + ((lsegkai/vect)**dimf(is, js)) ! incluye el jacobiano R(segmento)
           enddo
           enddo
         endif
         endif

      enddo!iz
      enddo!iy
      enddo!ix
      
      do jR = 1, dimR ! save all
      do jZ = -Xulimit, Xulimit
         do is=1,Npoorsv
         do js=1,Npoorsv
           Xu(ii, jR, jZ, is, js) = Xu(ii, jR, jZ, is, js)/(MCsteps**3)*(2.0*cutoff)**3
         enddo
         enddo
      enddo
      enddo
  end do ! ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end select 


endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do is=1,Npoorsv
do js=1,Npoorsv
  
  
  write(kaisfilename,'(A5,BZ,I3.3,A1,I3.3,A3)')'kais.',is,'.',js,'.in'
  
  open(unit=200,file='suma.dat')
  open(unit=is*110+js, file=kaisfilename)

  if (flagkai.eq.0) then

    read(is*110+js,*)nada
    read(is*110+js,*)curvkais,dimRkais,dimRinikais,Xulimitkais,dimfkais(is,js),MCfactorkais
  
    if (curvkais.ne.curvature) then
      print*,"curvature of kais non equal curvature of DEFINITIONS.txt"
      stop
    endif

    if (dimRkais.ne.dimR) then
      print*,"box size of kais non equal box size of DEFINITIONS.txt"
      stop
    endif

    if (Xulimitkais.ne.Xulimit) then
      print*,"Xulimit of kais non equal Xulimit of DEFINITIONS.txt"
      stop
    endif

    if (dimfkais(is,js).ne.dimf(is,js)) then
      print*,"dimf of kais non equal dimf of DEFINITIONS.txt"
      stop
    endif

    if (MCfactorkais.ne.MCfactor) then
      print*,"MCsteps prefactor of kais.XXX.XXX.in non equal to MCsteps prefactor defined in kai.f90"
      stop
    endif

    if (dimRinikais.ne.dimRini) then
      print*,"dimRini of kais non equal to dimRini of DEFINITIONS.txt"
      stop
    endif

  endif

  if (flagkai.eq.1) then

    write(is*110+js,*)'#curvature dimR dimRini Xulimit dimf factorMCsteps#'
    write(is*110+js,*)curvature,dimR,dimRini,Xulimit,dimf(is,js),MCfactor

  endif

  do ii=1, dimR
  do jR=1, dimR
  do jZ=-Xulimit,Xulimit

     if (flagkai.eq.1) then
       write(is*110+js,*)ii,jR,jZ,Xu(ii,jR,jZ,is,js) ! residual size of iteration vector
     endif

     if (flagkai.eq.0) then
  
       read(is*110+js,*)a,b,c,Xu(ii,jR,jZ,is,js)
  
       if (a.ne.ii) then
         print*,'a non equal ii'
         stop
       endif
 
       if (b.ne.jR) then
         print*,'b non equal jR'
         stop
       endif

       if (c.ne.jZ) then
         print*, 'c non equal jZ'
         stop
       endif
       
     endif

  enddo
  enddo
  enddo

  do jR = 1, dimR
  do jZ = -Xulimit,Xulimit
    sumaXu(is,js) = sumaXu(is,js) + Xu(1,jR,jZ,is,js)
  enddo
  enddo

  write(200,*)is,js,sumaXu(is,js)

  close(is*110+js)

  if (flagonekais.eq.1) then
    exit
  endif
enddo !js
if (flagonekais.eq.1) then
  exit
endif

enddo !is

if (flagonekais.eq.1) then
  do is=1,Npoorsv
  do js=1,Npoorsv
    sumaXu(is,js) = sumaXu(1,1)
    Xu(:,:,:,is,js) = Xu(:,:,:,1,1)
  enddo
  enddo
endif

close(200)
end

