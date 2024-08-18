subroutine creador

use globals
use layer
use volume
use seed1
use longs
use MPI
use transgauche
use mkai
use cadenaMD
implicit none

integer is

integer kk, jj,i,j,k,ii ! dummy indice0s

INTEGER temp_R
real*8 tempr_R, tempr_Z

integer ncha
!character*20 filename
real*8 chains(3,maxlong,ncha_max) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
real*8 Uconf
integer*1 Ntconf(maxlong)
real*8 sumUgyr, sumRgyr(0:Npoorsv+1), Rgyr(0:Npoorsv+1), Ugyr, Rgyrprom(0:Npoorsv+1)
!real*8 rog
integer conf              ! counts number of conformations


character*80 line
integer, external :: PBCSYMI


! MPI
integer tag
parameter(tag = 0)
integer NC
character*9 filename2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DEBUG
!do i = 1, long(1)
!write(filename,'(A4,I2.2,A4)')'seg.',i,'.dat'
!open(unit=9500+i,file=filename)
!write(filename,'(A4,I2.2,A4)')'sxg.',i,'.dat'
!open(unit=9600+i,file=filename)
!write(filename,'(A4,I2.2,A4)')'syg.',i,'.dat'
!open(unit=9700+i,file=filename)
!write(filename,'(A4,I2.2,A4)')'szg.',i,'.dat'
!open(unit=9800+i,file=filename)
!enddo

!open(unit=9900,file='rog.dat')
!open(unit=9917,file='e17.dat')
!open(unit=9957,file='e57.dat')
!open(unit=9915,file='e15.dat')
!open(unit=9927,file='e27.dat')
!open(unit=9947,file='e47.dat')




call initcha              ! init matrices for chain generation

innR = 0
innZ = 0
Uchain = 0.
Ntrans = 0

do NC = 1, Ncomp

sumRgyr=0.
sumUgyr=0.
Rgyrprom=0.
Uconf=0.
Ntconf=0.

conf=0                    ! counter of number of conformations

lineposMD = 0 

if (flagMD(NC).eq.1) then
    write(filename2,'(A3,I3.3,A3)')'MD.',NC,'.in'
    open(7777,file=filename2)
    if (rank.gt.0) then
       call MPI_RECV(lineposMD, 1, MPI_INTEGER, rank-1, rank-1, MPI_COMM_WORLD, status, ierr)
       print*, "rank", rank, "received from rank", rank-1,"linepos =", lineposMD

! Discard first linepos lines
       do i = 1, lineposMD
        read (7777, '(A)') line
       enddo

       call MPI_RECV(seed, 1, MPI_INTEGER, rank-1, rank-1, MPI_COMM_WORLD, status, ierr)
       print*, "rank", rank, "received from rank", rank-1,"seed =", seed
    endif
else
    if (rank.gt.0) then
       call MPI_RECV(seed, 1, MPI_INTEGER, rank-1, rank-1, MPI_COMM_WORLD, status, ierr)
       print*, "rank", rank, "received from rank", rank-1,"seed =", seed
    endif
endif


do while (conf.lt.cuantas(NC))

if(flagMD(NC).eq.1) then
   call cadenasMD(chains,ncha,Uconf,Ntconf,Ugyr,Rgyr,NC)
else 
   call cadenas(chains,ncha,Uconf,Ntconf,Ugyr,Rgyr,NC)
endif

!   open(unit=100000,file='chains.dat')
!   do i=1,3
!   do j=1,long(1)
!   do k=1,ncha
!     write(100000,*)conf+k-1,i,j,k,chains(i,j,k)
!   enddo
!   enddo
!   enddo


   do is=0,Npoorsv+1
      sumRgyr(is)=sumRgyr(is)+Rgyr(is)*exp(-Ugyr)
   enddo

   sumUgyr=sumUgyr+exp(-Ugyr)

   do j=1,ncha
      if(conf.lt.cuantas(NC)) then



! DEBUG
!      rog = 0.0
!      do k = 1, long(NC)
!      do kk = 1,long(NC)
!      rog = rog +(chains(1,k,j)-chains(1,kk,j))**2   &
!                +(chains(2,k,j)-chains(2,kk,j))**2   &
!                +(chains(3,k,j)-chains(3,kk,j))**2
!      enddo
!      enddo
!      rog = rog / float(long(NC)**2)/2.

!      write(9900,*)rog
!      write(9917,*)sqrt((chains(1,1,j)-chains(1,7,j))**2   &
!                +(chains(2,1,j)-chains(2,7,j))**2   &
!                +(chains(3,1,j)-chains(3,7,j))**2)

!      write(9957,*)sqrt((chains(1,5,j)-chains(1,7,j))**2   &
!                +(chains(2,5,j)-chains(2,7,j))**2   &
!                +(chains(3,5,j)-chains(3,7,j))**2)

!      write(9915,*)sqrt((chains(1,1,j)-chains(1,5,j))**2   &
!                +(chains(2,1,j)-chains(2,5,j))**2   &
!                +(chains(3,1,j)-chains(3,5,j))**2)

!      write(9927,*)sqrt((chains(1,2,j)-chains(1,7,j))**2   &
!                +(chains(2,2,j)-chains(2,7,j))**2   &
!                +(chains(3,2,j)-chains(3,7,j))**2)

!      write(9947,*)sqrt((chains(1,4,j)-chains(1,7,j))**2   &
!                +(chains(2,4,j)-chains(2,7,j))**2   &
!                +(chains(3,4,j)-chains(3,7,j))**2)

         conf=conf+1
         Uchain(conf,NC)=Uconf
         Ntrans(:,conf,NC) = Ntconf(:)
         do k=1,long(NC)
! DEBUG
!          write(9500+k,*)sqrt(chains(1,k,j)**2+chains(2,k,j)**2+chains(3,k,j)**2)
!          write(9600+k,*)chains(1,k,j)
!          write(9700+k,*)chains(2,k,j)
!          write(9800+k,*)chains(3,k,j)


            do ii = minntotR(NC),maxntotR(NC) ! position of first segment (or Center of mass?)

               select case (abs(curvature))
                 case (2)
                  tempr_R=((chains(1,k,j)+(float(ii+dimRini)-0.5)*deltaR)**2 + chains(2,k,j)**2 + chains(3,k,j)**2 )**(0.5)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case (1)
                  tempr_R=((chains(1,k,j)+(float(ii+dimRini)-0.5)*deltaR)**2+chains(2,k,j)**2)**(0.5)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case (0) 
                  tempr_R=abs(chains(1,k,j)+(float(ii+dimRini)-0.5)*deltaR)
                  temp_R=int(tempr_R/deltaR)+1  ! put them into the correct layer
                 case(3) ! lamella with PBC in r
                  tempr_R=chains(1,k,j)+(float(ii+dimRini)-0.5)*deltaR
                  temp_R=int(anint(tempr_R/deltaR)) 
                  temp_R = PBCSYMI(temp_R,dimR)+dimRini ! puts the segment within the calculation box using PBC,
                                                        ! adds dimRini for compatibility 
             endselect
             
               if(temp_R.gt.(dimR+dimRini)) then
                     print*, 'creador.f90: increase dimR', temp_R, chains(1,k,j), ii, j
                     do jj = 1, long(NC)
                        print*, jj, chains(:,jj,j)
                     enddo
                     stop
               endif
              innR(k,conf,ii,NC)=temp_R-dimRini ! in which layer is the segment "k" of a chain at position "ii" and conformation "conf"
            enddo ! ii
         tempr_Z=chains(3,k,j)
         innZ(k,conf,NC)=-int(anint(tempr_Z/deltaZ)) ! mirror conformation in Z (legacy, does not introduce a bias)
         enddo ! k
      endif
   enddo ! j

enddo ! while

if(flagMD(NC).eq.1) then
  close(7777)
  if (rank.lt.size-1) then
     call MPI_SEND(lineposMD, 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD, ierr)
     print*, "rank", rank, "sent to rank", rank+1,"linepos=", lineposMD

     call MPI_SEND(seed, 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD, ierr)
     print*, "rank", rank, "sent to rank", rank+1,"seed=", seed

  endif
else
  if (rank.lt.size-1) then
     call MPI_SEND(seed, 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD, ierr)
     print*, "rank", rank, "sent to rank", rank+1,"seed=", seed
  endif
endif

do is=0,Npoorsv+1
   Rgyrprom(is)=sumRgyr(is)/sumUgyr
enddo


call MPI_BARRIER(MPI_COMM_WORLD, ierr)


if(rank.eq.0) then
     print*," chains component ", NC, " out of ", Ncomp, " ready"
     do is=0,Npoorsv+1
       print*,is, Rgyrprom(is), sumRgyr(is), sumUgyr
     enddo
endif


! print Rgyr 

do is=0,Npoorsv+1
   write(2533+NC,*) is, Rgyrprom(is)
enddo
close(2533+NC)

enddo ! NC

! calc vchain, moved here because need segpoorsv

vchain=0.0
do NC = 1,Ncomp
do i=1,long(NC)
  vchain(NC)=vchain(NC)+vpol(segpoorsv(i,NC))
enddo
enddo

!stop
end
