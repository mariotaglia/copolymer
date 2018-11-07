subroutine creador

use globals
use layer
use volume
use seed1
use longs
use MPI
use transgauche

implicit none

real*8 solvetime1, solvetime2, solveduration
integer is,ic
integer *4 ier ! Kinsol error flag
real*8 pi
real*8 Na               
parameter (Na=6.02d23)
integer av1(ntot), av2(ntot)
real*8 avtmp

integer i,j,k,m,ii,flag,c, jj, iR,iZ ! dummy indice0s

INTEGER temp_R, temp_Z
real*8 tempr_R, tempr_Z
real*8 tmp

real*8 min1               ! variable to determine minimal position of chain
integer qqq,www,eee

integer il,inda,ncha

REAL*8 xfile((npoorsv+2)*ntot)                        
real*8 algo, algo2                  

real*8 chains(3,long,ncha_max) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
real*8 zp(long)
real*8 Uconf
integer*1 Ntconf(long)
real*8 sum,sumel          ! auxiliary variable used in free energy computation  
real*8 sumpi,sumrho,sumrhopol, sumrho2, sumrho2mol !suma de la fraccion de polimero
real*8 sumUgyr, sumRgyr(0:Npoorsv+1), Rgyr(0:Npoorsv+1), Ugyr, Rgyrprom(0:Npoorsv+1)

integer conf              ! counts number of conformations
INTEGER cc

! MPI
integer tag, source
parameter(tag = 0)
integer err
integer ier_tosend
double  precision norma_tosend

integer in1tmp(long,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call initcha              ! init matrices for chain generation

conf=0                    ! counter of number of conformations

if (rank.gt.0) then
   call MPI_RECV(seed, 1, MPI_INTEGER, rank-1, rank-1, MPI_COMM_WORLD, status, ierr)
   print*, "rank", rank, "received from rank", rank-1,"seed =", seed
endif

innR = 0
innZ = 0
sumRgyr(:)=0.
sumUgyr=0.
Rgyrprom(:)=0.

do while (conf.lt.cuantas)

   call cadenas(chains,ncha,Uconf,Ntconf,Ugyr,Rgyr)
   
   do is=0,Npoorsv+1
      sumRgyr(is)=sumRgyr(is)+Rgyr(is)*exp(-Ugyr)
   enddo

   sumUgyr=sumUgyr+exp(-Ugyr)

   do j=1,ncha

      if(conf.lt.cuantas) then

         conf=conf+1
         Uchain(conf)=Uconf
         Ntrans(:,conf) = Ntconf(:)

         do k=1,long
            do ii = 1,maxntotR ! position of first segment (or Center of mass?)

               select case (abs(curvature))
                 case (2)
                  tempr_R=((chains(1,k,j)+(float(ii)-0.5)*deltaR)**2 + chains(2,k,j)**2 + chains(3,k,j)**2 )**(0.5)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case (1)
                  tempr_R=((chains(1,k,j)+(float(ii)-0.5)*deltaR)**2+chains(2,k,j)**2)**(0.5)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case (0) 
                  tempr_R=abs(chains(1,k,j)+(float(ii)-0.5)*deltaR)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
               endselect
             
               if(temp_R.gt.dimR) then
                  if(rank.eq.0)print*, 'main.f90: increase dimR'
                  stop
               endif

              innR(k,conf,ii)=temp_R ! in which layer is the segment "k" of a chain at position "ii" and conformation "conf"
         
            enddo ! ii
         
         tempr_Z=chains(3,k,j)
         innZ(k,conf)=int(anint(tempr_Z/deltaZ))

         enddo ! k

      endif

   enddo ! j

enddo ! while

if (rank.lt.size-1) then
   call MPI_SEND(seed, 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD, ierr)
   print*, "rank", rank, "sent to rank", rank+1,"seed=", seed
endif

do is=0,Npoorsv+1
   Rgyrprom(is)=sumRgyr(is)/sumUgyr
enddo

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if(rank.eq.0) then

   print*," chains ready"

   do is=0,Npoorsv+1
      print*,is, Rgyrprom(is), sumRgyr(is), sumUgyr
   enddo

endif
end
