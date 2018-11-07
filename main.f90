

call initmpi
call parser
call init
call allocation
call kai
call creador
call solve

end

subroutine solve

use globals
use MPI
use mkai
implicit none

integer is
integer *4 ier ! Kinsol error flag
real*8 x1((npoorsv+2)*ntot),xg1((npoorsv+2)*ntot),x1ini((npoorsv+2)*ntot)   ! density solvent iteration vector

real*8 trash, trash2

integer n                 ! number of lattice sites

external fcnelect         ! function containing the SCMFT eqs for solver
integer i,iR,iZ ! dummy indices

REAL*8 xfile((npoorsv+2)*ntot)                        

integer countfile         ! enumerates the outputfiles 


! MPI
integer tag
parameter(tag = 0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     initial guess all 1.0 

n=dimR*dimZ

do iR=1,dimR
do iZ=1,dimZ

   xg1(dimR*(iZ-1)+iR)=1.0
   x1(dimR*(iZ-1)+iR)=1.0

   do is=1,Npoorsv+1
      xg1(n*is+dimR*(iZ-1)+iR)=0.0001
      x1(n*is+dimR*(iZ-1)+iR)=0.0001
   enddo

enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     init guess from files fort.100 (solvent) and fort.200 (potential)                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (infile.ge.1) then
   do iR=1,dimR
   do iZ=1,dimZ

     read(100,*)trash,trash2,xfile(dimR*(iZ-1)+iR)   ! solvent
     x1(dimR*(iZ-1)+iR)=xfile(dimR*(iZ-1)+iR)
     xg1(dimR*(iZ-1)+iR)=xfile(dimR*(iZ-1)+iR)
     do is=1,Npoorsv 

       read(100+is,*)trash,trash2,xfile(n*is+dimR*(iZ-1)+iR)   ! poorsolvent desde 1 a npoorsv 
!       if(xfile(i+n*is).lt.1.0d-30)xfile(i+n*is)=1.0d-30
       x1(n*is+dimR*(iZ-1)+iR)=xfile(n*is+dimR*(iZ-1)+iR)
       xg1(n*is+dimR*(iZ-1)+iR)=xfile(n*is+dimR*(iZ-1)+iR)

     enddo !is

     read(200,*)trash,trash2,xfile(n*(Npoorsv+1)+dimR*(iZ-1)+iR)
     x1(n*(npoorsv+1)+dimR*(iZ-1)+iR)=xfile(n*(npoorsv+1)+dimR*(iZ-1)+iR)
     xg1(n*(npoorsv+1)+dimR*(iZ-1)+iR)=xfile(n*(npoorsv+1)+dimR*(iZ-1)+iR)

   enddo !iR
   enddo !iZ

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     computation starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     initializations of input depended variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


iter=0                    ! iteration counter
countfile=1

actionflag = 0 ! Actionflag controls the current action of loop
               ! = 0 loop maxntotcounter from 2 to maxntot
               ! = 1 increases npol from npolini to npollast
               ! = 2 decreases npol from npolini to npolfirst
               ! = 3 finalize

maxntotcounterR = maxntotR
maxntotcounterZ = maxntotZ
npol = npolini

do while (actionflag.lt.3)
   if(rank.eq.0)print*, ' npol:', npol, 'maxntotR:', maxntotcounterR, 'maxntotZ:', maxntotcounterZ

   do i=1,(npoorsv+2)*n             ! initial guess for x1
      xg1(i)=x1(i)
   enddo

! JEFE
   if(rank.eq.0) then ! solo el jefe llama al solver
      print*, 'solve: Enter solver ', (npoorsv+2)*ntot, ' eqs'
   endif   

   iter=0
   call call_kinsol(x1, xg1, ier)

   do iR=1,dimR
   do iZ=1,dimZ
     xsol(iR,iZ)=x1(dimR*(iZ-1)+iR)
   enddo
   enddo

   if((norma.gt.error).or.(ier.lt.0).or.(isnan(norma))) then
      if(actionflag.gt.0) then
         print*, " I am ", rank, " I stopped the work"
         stop
      endif
   endif

   print*, "norma", rank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calc free energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call calc_free_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! save to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(rank.eq.0)call save2disk(actionflag, countfile)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  next calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  countfile = countfile+1 ! next

   select case (actionflag)

     case(0) ! maxntot loop
      if (maxntotcounterR.lt.maxntotR) then
        maxntotcounterR=maxntotcounterR+1
      endif

      if(maxntotcounterR.eq.maxntotR) then
        actionflag=1 
        countfile = 1
      endif

     case(1)  ! increases from npolini to npollast

      if(npol.eq.npolini)x1ini=x1
      npol = npol + npolstep      
      if(npol.gt.npollast) then
      npol = npolini
      actionflag = 2
      countfile = 1
      x1 = x1ini

      endif       

     case(2)
      npol = npol - npolstep
      if(npol.lt.npolfirst)actionflag = 3

   endselect

end do ! loop de npol

!!!!!!!!!!!!!!!!!!!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MPI_FINALIZE(ierr) ! finaliza MPI

stop

end
