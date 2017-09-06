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
integer seed
real*8 xmin,xmax,ymin,ymax,zmin,zmax
integer MCsteps ! numero de steps de MC


real*8 R,theta,z
real*8 rn
integer i, ii, is, js, a, b
real*8 rands
real*8 pi
real*8 x1,x2,y1, y2, z1, z2, vect
integer iR, ix,iy,iz, itheta
integer j
real*8 radio
real*8 cutoff
real*8, allocatable :: sumaXu(:,:)
character*16 kaisfilename

if(rank.eq.0)print*,'Kai calculation'

allocate(sumaXu(Npoorsv,Npoorsv))

cutoff = (float(Xulimit)+0.5)*delta

pi=dacos(-1.0d0)          ! pi = arccos(-1) 
radio = float(ntot)*delta

Xu = 0.0 ! vector Xu

seed = 1010
MCsteps = 60*Xulimit
sumaXu(:,:)=0.0

if (flagkai.eq.1) then

do ii = 1, ntot ! loop sobre cada posicion del segmento

      ymax = cutoff
      ymin = -cutoff

      zmax = cutoff
      zmin = -cutoff

      xmax = cutoff + (dfloat(ii) - 0.5)*delta
      xmin = -cutoff + (dfloat(ii) - 0.5)*delta
      
      do ix = 1, MCsteps
      do iy = 1, MCsteps
      do iz = 1, MCsteps

! coordenadas del segmento (x1,y1,z1) y del punto a integrar (x2,y2,z2)

         x1 = (dfloat(ii) - 0.5)*delta ! asume theta segmento = 0, z segmento = 0 y segmento en el centro de la layer
         y1 = 0.0
         z1 = 0.0

         x2 = xmin + (xmax-xmin)*dfloat(ix-1)/dfloat(MCsteps-1)
         y2 = ymin + (ymax-ymin)*dfloat(iy-1)/dfloat(MCsteps-1)
         z2 = zmin + (zmax-zmin)*dfloat(iz-1)/dfloat(MCsteps-1)

         select case (abs(curvature))
         case (0)
         R = abs(x2)
         case (1)
         R = sqrt(x2**2 + y2**2)
         case (2)
         R = sqrt(x2**2 + y2**2 + z2**2)
         end select

         vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diferencia
         j = int(R/delta)+1 ! j tiene la celda donde cae el punto a integrar

         if(j.le.ntot) then


         if(vect.le.(cutoff)) then ! esta dentro de la esfera del cut-off   
         if(vect.ge.lseg) then ! esta dentro de la esfera del segmento
           do is=1,Npoorsv
           do js=1,Npoorsv
              Xu(ii, j, is, js) = Xu(ii, j, is, js) + ((lseg/vect)**dimf(is, js)) ! incluye el jacobiano R(segmento)
           enddo
           enddo
         endif
         endif

         endif
         
      enddo!iz
      enddo!iy
      enddo!ix
      
      do j = 1, ntot
         do is=1,Npoorsv
         do js=1,Npoorsv
         Xu(ii, j,is,js) = Xu(ii, j, is, js)/(MCsteps**3)*(2.0*cutoff)**3
         enddo
         enddo
      enddo
end do ! ii

endif

do is=1,Npoorsv
do js=1,Npoorsv

write(kaisfilename,'(A5,BZ,I3.3,A1,I3.3,A4)')'kais.',is,'.',js,'.dat'
open(unit=200,file='suma.dat')
open(unit=is*110+js, file=kaisfilename)


  do ii=1,ntot
  do j=1,ntot

     if (flagkai.eq.1) then
     write(is*110+js,*)ii,j,Xu(ii,j,is,js) ! residual size of iteration vector
     endif

     if (flagkai.eq.0) then
     read(is*110+js,*)a,b,Xu(ii,j,is,js)

      if (a.ne.ii) then
      print*,'a non equal ii'
      stop
      endif
 
      if (b.ne.j) then
      print*,'b non equal j'
      stop
      endif

     endif

  enddo
  enddo

do i = 20-Xulimit, 20+Xulimit
sumaXu(is,js) = sumaXu(is,js) + Xu(20,i,is,js)
enddo

write(200,*)is,js,sumaXu(is,js)

close(is*110+js)

enddo !js
enddo !is

close(200)
end

