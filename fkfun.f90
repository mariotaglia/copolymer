subroutine fkfun(x,f,ier2)
use globals
use mcharge
use partfunc
use layer
use volume
use bulk
use longs
use MPI
use pis
use mkai
use transgauche
implicit none

real*8 all_tosend(4*ntot), all_toreceive(4*ntot)
integer*4 ier2
real*8 protemp, sttemp
real*8 x((Npoorsv+2)*ntot),f((Npoorsv+2)*ntot)
real*8 xh(dimR+1,dimZ) 
real*8 xpot(0:Npoorsv,dimR,dimZ), xpot_a(Nacids,dimR,dimZ), xpot_b(Nbasics,dimR,dimZ)
real*8 pro(cuantas)
real*8 time1, time2, duration, looptime1, looptime2, loopduration
integer iR,iZ,kZ,kkZ,k,i,j,k1,k2,ii,jj,ic,aR,aZ,iZm,iZp,jZp,jZm,as,bs        ! dummy indices
integer is, js,ia,ib,iiR,iiZ,jR,jZ
integer err
integer n
real*8 avpol_tmp(0:Npoorsv,dimR,dimZ), avpola_tmp(Nacids,dimR,dimZ), avpolb_tmp(Nbasics,dimR,dimZ)
real*8 avpol_tosend(0:Npoorsv,dimR,dimZ), avpola_tosend(Nacids,dimR,dimZ), avpolb_tosend(Nbasics,dimR,dimZ)
real*8 xpol_tosend(dimR,dimZ)
real*8 algo, algo1,algo2
double precision, external :: factorcurv
real*8 sumpol
real*8 q_tosend(dimR,dimZ)
real*8 sumprolnpro_tosend(dimR,dimZ), sumprouchain_tosend(dimR,dimZ)
real*8 sumtrans_tosend(dimR,dimZ,long)
real*8 sumtrans(dimR,dimZ,long)
real*8 gradphi2
integer, external :: PBCSYMI
! Jefe
!flagsolver=1

!if(rank.eq.0) then ! llama a subordinados y pasa vector x
!   time1=MPI_WTIME()
!   CALL MPI_BCAST(x, (Npoorsv+1)*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
!endif

n = ntot 

! Recover xh and phi from input

do iR=1,dimR
do iZ=1,dimZ
  xh(iR,iZ)=x(dimR*(iZ-1)+iR) ! solvent volume fraction (xh) ir read from the x provided by kinsol
enddo
enddo !OJO boundary conditions?

xh(iR+1,iZ)=1.0


do iR=1,dimR
do iZ=1,dimZ
  phi(iR,iZ)=x(n*(Npoorsv+1)+dimR*(iZ-1)+iR)
enddo
enddo

phi(0,:)=phi(1,:)
phi(dimR+1,:)=0.0


! Recover xtotal from input
do is = 1,Npoorsv
  do iR=1,dimR
  do iZ=1,dimZ
    xtotal(is,iR,iZ) = x(n*is+dimR*(iZ-1)+iR) !  segments volume fractions (xtotal) are read from the x provided by kinsol
  enddo
  enddo
enddo

avpos=0.0
avneg=0.0

do iR=1,dimR
do iZ=1,dimZ

   avpos(iR,iZ)=vpos*expmupos*xh(iR,iZ)**vpos*dexp(-phi(iR,iZ)) ! volume fraction of cations
   avneg(iR,iZ)=vneg*expmuneg*xh(iR,iZ)**vneg*dexp(phi(iR,iZ)) ! volume fraction of anions
   avHplus(iR,iZ)=expmuHplus*xh(iR,iZ)*dexp(-phi(iR,iZ)) ! volume fraction of H+
   avOHmin(iR,iZ)=expmuOHmin*xh(iR,iZ)*dexp(phi(iR,iZ)) ! volume fraction of OH-

   do is=1,Nacids
     fAmin(is,iR,iZ)=1.0/(avHplus(iR,iZ)/(Ka(is)*xh(iR,iZ))+1.0)
   enddo

   do is=1,Nbasics
     fBHplus(is,iR,iZ)=1.0/(avOHmin(iR,iZ)/(Kb(is)*xh(iR,iZ))+1.0)
   enddo

enddo
enddo


! Caculation of dielectric function
! Everything that it is not water or ions has dielectric dielP

do iR = 1, dimR
do iZ = 1, dimZ
  dielpol(iR,iZ) = 1.0 - xh(iR,iZ) - avpos(iR,iZ) - avneg(iR,iZ) - avHplus(iR,iZ) - avOHmin(iR,iZ)
enddo
enddo
call dielectfcn(dielpol,epsfcn,Depsfcn)



! Calculation of xpot

do iZ = 1, dimZ
do iR = 1, dimR
! osmotic pressure
   xpot(0,iR,iZ) = xh(iR,iZ)**vpol(0) ! exp(-pi(r)v_pol) / units of v_pol: nm^3
enddo   
enddo

! dielectrics
do iZ = 1,dimZ
   jZp=iZ+1
   iZp=PBCSYMI(jZ,dimZ)
do iR = 1,dimR
   gradphi2 = ((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2+((phi(iR,iZp)-phi(iR,iZ))/deltaZ)**2
   xpot(0,iR,iZ) = xpot(0,iR,iZ)*exp(Depsfcn(iR,iZ)*gradphi2*vpol(0)*vsol*wperm/2.0)
enddo 
enddo

do iZ = 1, dimZ
do iR = 1, dimR

  do is = 1, Npoorsv
!   calculate xpot(i, is)

    protemp = 0.0
    
      do js = 1, Npoorsv 
         do jR = 1, dimR
         do jZ = -Xulimit, Xulimit
            kZ=jZ+iZ
            kkZ=PBCSYMI(kZ,dimZ)
            protemp = protemp+st(is,js)*Xu(iR,jR,jZ,is,js)*xtotal(js,jR,kkZ)/(vpol(js)*vsol) ! vpol*vsol in units of nm^3
         enddo
         enddo
      enddo

    xpot(is,iR,iZ) = xh(iR,iZ)**vpol(is) ! exp(-pi(r)v_pol) / units of v_pol: nm^3
    xpot(is,iR,iZ) = xpot(is,iR,iZ)*dexp(protemp) 

  enddo

enddo
enddo

do is= 1, Npoorsv
  do iZ= 1, dimZ
  jZp=iZ+1
  iZp=PBCSYMI(jZp,dimZ)
    do iR= 1, dimR
      gradphi2 = ((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2+((phi(iR,iZp)-phi(iR,iZ))/deltaZ)**2
      xpot(is,iR,iZ) = xpot(is,iR,iZ)*exp(Depsfcn(iR,iZ)*gradphi2*vpol(is)*vsol*wperm/2.0)
    enddo
  enddo
enddo

do iR = 1, dimR
do iZ = 1, dimZ
  do ic = 1,Nacids
    xpot_a(ic,iR,iZ) = 1.0/fAmin(ic,iR,iZ)*exp(phi(iR,iZ))
  enddo
  do ic = 1,Nbasics
    xpot_b(ic,iR,iZ)= 1.0/fBHplus(ic,iR,iZ)*exp(-phi(iR,iZ))
  enddo
enddo
enddo

!    probability distribution

avpola_tosend = 0.0
avpola_tmp = 0.0
avpola = 0.0
avpolb_tosend = 0.0
avpolb_tmp = 0.0
avpolb = 0.0
avpol_tosend = 0.0
xpol_tosend = 0.0
avpol_tmp = 0.0
avpol = 0.0
xpol = 0.0
q = 0.0
q_tosend=0.0d0                   ! init q to zero
sumprolnpro_tosend = 0.0
sumprolnpro = 0.0
sumprouchain_tosend=0.0
sumprouchain=0.0
sumtrans_tosend = 0.0
sumtrans = 0.0
all_tosend = 0.0
all_toreceive = 0.0


do iiR=1, maxntotcounterR ! position of center of mass 
do iiZ=1, maxntotcounterZ
   do i=1, cuantas ! loop over conformations
 
      pro(i) = exp(-Uchain(i))

       do k=1,long
         
         iZ = innZ(k,i)+iiZ
         aZ = PBCSYMI(iZ,dimZ)
         aR = innR(k,i,iiR)
         is = segpoorsv(k)
         ia = acidtype(k)
         ib = basictype(k) 

         pro(i)= pro(i) * xpot(is,aR,aZ)

         pro(i)= pro(i) * xpot_a(ic,aR,aZ) 
     
         pro(i)= pro(i) * xpot_b(ic,aR,aZ)
 
      enddo !k

      q_tosend(iiR,iiZ)=q_tosend(iiR,iiZ)+pro(i) ! all_tosend(ii) = all_tosend(ii) + pro(i) 
      sumprolnpro_tosend(iiR,iiZ) = sumprolnpro_tosend(iiR,iiZ) + pro(i)*dlog(pro(i)) ! all_tosend(ntot+ii) = all_tosend(ntot+ii) + pro(i)*dlog(pro(i))
      sumprouchain_tosend(iiR,iiZ) = sumprouchain_tosend(iiR,iiZ) + pro(i)*Uchain(i) !all_tosend(ntot*2+ii) = all_tosend(ntot*2+ii) + pro(i)*Uchain(i) 
      xpol_tosend(iiR,iiZ) = xpol_tosend(iiR,iiZ)+pro(i) ! all_tosend(ntot*3+ii) = all_tosend(ntot*3+ii) + pro(i)


      do j = 1, long ! loop over number of segments

         iZ = innZ(k,i)+iiZ
         aZ = PBCSYMI(iZ,dimZ)
         aR = innR(k,i,iiR)
         is = segpoorsv (j)
         as = acidtype (j)
         bs = basictype (j)

         sumtrans_tosend(iiR,iiZ,j) =  sumtrans_tosend(iiR,iiZ,j) +  pro(i)*float(Ntrans(j,i))
         avpol_tmp(is,aR,aZ)=avpol_tmp(is,aR,aZ)+pro(i)*factorcurv(iiR,aR) ! avpol_tmp is avg number of segments "is" at position "j" 
         avpola_tmp(ia,aR,aZ)=avpola_tmp(ia,aR,aZ)+pro(i)*factorcurv(iiR,aR) ! avpola_tmp is avg number of acid segments "ic" at position "j"
         avpolb_tmp(ib,aR,aZ)=avpolb_tmp(ib,aR,aZ)+pro(i)*factorcurv(iiR,aR) ! avpolb_tmp is avg number of basic segments "ic" at position "j" 

      enddo ! j

   enddo ! i

enddo ! iiR
enddo ! iiZ

avpol_tosend(:, 1:dimR, 1:dimZ)=avpol_tmp(:, 1:dimR, 1:dimZ) 
avpola_tosend(:, 1:dimR, 1:dimZ)=avpola_tmp(:, 1:dimR, 1:dimZ)
avpolb_tosend(:, 1:dimR, 1:dimZ)=avpolb_tmp(:, 1:dimR, 1:dimZ)

!------------------ MPI -----------------`-----------------------------
!1. Todos al jefe

!call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe

!  Junta avpol       

! Subordinados


!  Junta avpol       

   call MPI_ALLREDUCE(avpola_tosend, avpola, Nacids*ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(avpolb_tosend, avpolb, Nbasics*ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(avpol_tosend, avpol, (Npoorsv+1)*ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(xpol_tosend, xpol, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
!   call MPI_ALLREDUCE(all_tosend, all_toreceive, ntot*4, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(q_tosend, q, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(sumprolnpro_tosend, sumprolnpro, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(sumprouchain_tosend, sumprouchain, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(sumtrans_tosend, sumtrans, ntot*long, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!if(rank.ne.0) then
!   ier2 = 0
!   goto 3333

!endif

!if (rank.eq.0) then
!   looptime2=MPI_WTIME()
!   loopduration=looptime2-looptime1
!endif

! norma avpol
! integrate over whole system

! q(:) = all_toreceive(1:ntot)
! sumprolnpro(:) = all_toreceive(ntot+1:2*ntot)
! sumprouchain(:) = all_toreceive(ntot*2+1:3*ntot)
! xpol(:) = all_toreceive(ntot*3+1:4*ntot)


sumpol = 0.0

do iR = 1, dimZ
do iZ = 1, dimR
   do is = 0, Npoorsv

      select case (curvature)

       case (0)
        sumpol = sumpol + avpol(is,iR,iZ)*deltaR*deltaZ ! final result in units of chains/nm^2 (1D) or in units of chains/nm of belt (2D)
       case(1)
        sumpol = sumpol + avpol(is,iR,iZ)*deltaZ*(float(iR)-0.5)*deltaR*deltaR*2.0*pi ! final result in units of chains/nm of fiber (1D) or chains/fiber (2D)
       case(2)
        sumpol = sumpol + avpol(is,iR,iZ)*(((float(iR)-0.5)*deltaR)**2)*deltaR*4.0*pi ! final result in units of chains/micelle

      end select

   enddo
enddo
enddo

sumpol = sumpol/(vchain*vsol) 
avpol = avpol/sumpol*npol ! integral of avpol is fixed
avpola = avpola/sumpol*npol
avpolb = avpolb/sumpol*npol
sumpol = 0.0


do iR = 1, dimR
do iZ = 1, dimZ
   select case (curvature)
    case (0)
     sumpol = sumpol + xpol(iR,iZ)*deltaR*deltaZ ! final result in units of chains/nm^2 (1D) or chains/nm of belt (2D)
    case (1)
     sumpol = sumpol + xpol(iR,iZ)*deltaZ*(float(iR)-0.5)*deltaR*deltaR*2.0*pi ! final result in units of chains/nm (1D) or chains/fiber (2D)
    case (2)
     sumpol = sumpol + xpol(iR,iZ)*(((float(iR)-0.5)*deltaR)**2)*deltaR*4.0*pi ! final result in units of chains/micelle
   end select
enddo
enddo

xpol = xpol/sumpol*npol ! integral of avpol is fixed

trans = 0.0

do iR = 1, maxntotcounterR
do iZ = 1, maxntotcounterZ
   select case (curvature)
    case (0)
     trans(:) = trans(:) + sumtrans(iR,iZ,:)/q(iR,iZ)*xpol(iR,iZ)*deltaR*deltaZ ! final result in units of chains/nm^2 (1D) or chains/nm of belt (2D)
    case(1)
     trans(:) = trans(:) + sumtrans(iR,iZ,:)/q(iR,iZ)*xpol(iR,iZ)*deltaZ*(float(iR)-0.5)*deltaR*deltaR*2.0*pi ! final result in units of chains/nm (1D) or chains/fiber (2D)
    case(2)
     trans(:) = trans(:) + sumtrans(iR,iZ,:)/q(iR,iZ)*xpol(iR,iZ)*(((float(iR)-0.5)*deltaR)**2)*deltaR*4.0*pi ! final result in units of chains/micelle
   end select
enddo
enddo

trans(:) = trans(:)/npol

! contruction of f and the volume fractions

do iR = 1, dimR
do iZ = 1, dimZ

   f(dimR*(iZ-1)+iR)=xh(iR,iZ)+avneg(iR,iZ)+avpos(iR,iZ)+avHplus(iR,iZ)+avOHmin(iR,iZ)-1.0d0

   do is=0, Npoorsv
      f(dimR*(iZ-1)+iR) = f(dimR*(iZ-1)+iR) + avpol(is,iR,iZ)
   enddo

enddo
enddo

xcharge(:,:) = avpos(:,:)/(vpos*vsol)-avneg(:,:)/(vneg*vsol)+avHplus(:,:)/vsol-avOHmin(:,:)/vsol ! xcharge is avg charge density


do iR = 1, dimR
do iZ = 1, dimZ

   do ic= 1,Nacids
     xcharge(iR,iZ)=xcharge(iR,iZ)-avpola(ic,iR,iZ)*fAmin(ic,iR,iZ)/(vpol_a(ic)*vsol)
   enddo

   do ic= 1,Nbasics
     xcharge(iR,iZ)=xcharge(iR,iZ)+avpolb(ic,iR,iZ)*fBHplus(ic,iR,iZ)/(vpol_b(ic)*vsol)
   enddo
   iZp=iZ+1
   iZm=iZ-1
   jZp = PBCSYMI(iZp,dimZ)
   jZm = PBCSYMI(iZm,dimZ)
   select case (curvature)

    case (0)
     f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=xcharge(iR,iZ) &
     +wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-2.0*phi(iR,iZ)+phi(iR-1,iZ))*deltaR**(-2) &
     +wperm*epsfcn(iR,iZ)*(phi(iR,jZp)-2.0*phi(iR,iZ)+phi(iR,jZm))*deltaZ**(-2) &
     +wperm*(epsfcn(iR+1,iZ)-epsfcn(iR,iZ))*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2) &
     +wperm*(epsfcn(iR,jZp)-epsfcn(iR,iZ))*(phi(iR,jZp)-phi(iR,iZ))*deltaZ**(-2) 

    case(1)
     f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=xcharge(iR,iZ) &
     + wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2)/(float(iR)-0.5) &
     +wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-2.0*phi(iR,iZ)+phi(iR-1,iZ))*deltaR**(-2) &
     +wperm*epsfcn(iR,iZ)*(phi(iR,jZp)-2.0*phi(iR,iZ)+phi(iR,jZm))*deltaZ**(-2) &
     +wperm*(epsfcn(iR+1,iZ)-epsfcn(iR,iZ))*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2) &
     +wperm*(epsfcn(iR,jZp)-epsfcn(iR,iZ))*(phi(iR,jZp)-phi(iR,iZ))*deltaZ**(-2)

    case(2)
     f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=xcharge(iR,iZ) &
     + 2.0*wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2)/(float(iR)-0.5) &
     +wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-2.0*phi(iR,iZ)+phi(iR-1,iZ))*deltaR**(-2) &
     +wperm*(epsfcn(iR+1,iZ)-epsfcn(iR,iZ))*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2)

    end select

   f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)/(-2.0)

enddo
enddo

do is=1,Npoorsv
   do iR=1,dimR
   do iZ=1,dimZ ! xtotal
      f(n*is+dimR*(iZ-1)+iR) = -avpol(is,iR,iZ)+xtotal(is,iR,iZ) 
   enddo
   enddo
enddo


iter=iter+1

algo = 0.0
algo1 = 0.0
algo2 = 0.0

do i = 1, n*(Npoorsv+1)
   algo = algo + f(i)**2
end do

if(rank.eq.0) then 
   print*, iter, algo, q(1,1), Q(2,1), q(3,1), q(4,1)
endif

norma=algo

!3333 continue


ier2 = 0

return
end

! OJO 
! NOTAS MARIO SOBRE EQUIVALENCIA ENTRE WPERM Y CONSTQ
!
! 
!segun xpot:

!1/constq*vpol = wperm*vpol*vsol/delta**2
![wperm] = e^2/(kBT.nm) 
!psi en unidades de kBT/(e)
!constq = delta**2 / vsol / wperm

!segun fs:

!xcharge*vsol = qtot
!xcharge/wperm * delta^2 = qtot*constq
!1/wperm = vol*constq/delta^2

!segun input:

!constq = delta^2 * 4 pi / vsol * lb 
!constq = delta^2 * / vsol / wperm 

!segun free-energy

!1/constq = wperm/delta**2 * vsol
!constq = delta**2 / vsol / wperm

