!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 1D-coordenadas
! polares usando un metodo de MC
!
!#####################################################################


subroutine kai

implicit none
include global
include layer
include volume 

real*8 suma(ntot, ntot)
integer seed

integer MCsteps ! numero de steps de MC

real*8 R,theta,z
real*8 rn
integer i, ii
real*8 rands
real*8 pi
real*8 x1,x2,y1, y2, z1, z2, vect
integer iR, iz, itheta
integer j
real*8 radio

print*,'Kai calculation'
open(unit=111, file='kais.dat')

pi=dacos(-1.0d0)          ! pi = arccos(-1) 
radio = float(ntot)*delta

suma = 0.0

Xu = 0.0 ! vector Xu

seed = 1010
MCsteps = 20
!      MCsteps = 1

do ii = 1, ntot ! loop sobre cada posicion del segmento

do itheta = 1, MCsteps
do iz = 1, MCsteps
do iR = 1, MCsteps

theta = 2*pi*dfloat(itheta - 1)/dfloat(MCsteps)    ! numero 0 y 2pi
z = 3.0*(dfloat(iz - 1)/dfloat(MCsteps)-0.5)*delta ! numero entre -1.5*delta y 1.5*delta
R = radio*(dfloat(iR-1)/dfloat(MCsteps))         ! numero 0 y radio

! coordenadas del segmento (x1,y1,z1) y del punto a integrar (x2,y2,z2)

x1 = (dfloat(ii) - 0.5)*delta ! asume theta segmento = 0, z segmento = 0 y segmento en el centro de la layer
y1 = 0.0
z1 = 0.0
x2 = R*cos(theta)      y2 = R*sin(theta)
z2 = z

vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diferencia
j = int(R/delta)+1 ! j tiene la celda donde cae el punto a integrar

suma(ii, j) = suma(ii, j) + R

if(vect.le.(1.5*delta)) then ! esta dentro de la esfera del cut-off   
if(vect.ge.lseg) then ! esta dentro de la esfera del segmento
Xu(ii, j) = Xu(ii, j) + ((lseg/vect)**6)*R ! incluye el jacobiano R(segmento)
endif
endif

enddo
enddo
enddo


do j = 1, ntot
Xu(ii, j) = Xu(ii, j)/MCsteps**3*(3.0*delta)*2*pi*radio
suma(ii, j) = suma(ii, j)/MCsteps**3*(3.0*delta)*2*pi*(radio)
write(111,*)ii,j,Xu(ii,j) ! residual size of iteration vector
enddo

end do ! ii

close(111)

end

