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
real*8 xh(2*ntot) 
real*8 xpot(0:Npoorsv,2*ntot), xpot_a(Nacids,2*ntot), xpot_b(Nbasics,2*ntot)
real*8 pro(cuantas)
real*8 time1, time2, duration, looptime1, looptime2, loopduration
integer k,i,j,k1,k2,ii, jj,iz,ic       ! dummy indices
integer is, js
integer err
integer n
real*8 avpol_tmp(0:Npoorsv,2*ntot), avpola_tmp(Nacids,2*ntot), avpolb_tmp(Nbasics,2*ntot)
real*8 avpol_tosend(0:Npoorsv, ntot), avpola_tosend(Nacids,ntot), avpolb_tosend(Nbasics,ntot)
real*8 xpol_tosend(ntot)
real*8 algo, algo1,algo2
double precision, external :: factorcurv
real*8 sumpol
real*8 q_tosend(ntot)
real*8 sumprolnpro_tosend(ntot), sumprouchain_tosend(ntot)
real*8 sumtrans_tosend(ntot,long)
real*8 sumtrans(ntot,long)
real*8 inverse_of_vpolvsol, inverse_of_two
real*8 gradphi2

! Jefe
!flagsolver=1

!if(rank.eq.0) then ! llama a subordinados y pasa vector x
!   time1=MPI_WTIME()
!   CALL MPI_BCAST(x, (Npoorsv+1)*ntot , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
!endif

n = ntot

! Recover xh and phi from input
xh(1:n)=x(1:n) ! solvent volume fraction (xh) ir read from the x provided by kinsol
xh(n+1:2*n)=xsolbulk

phi(1:n)=x(n*(Npoorsv+1)+1:(Npoorsv+2)*n) ! electric potential (phi) is read from the x provided by kinsol
phi(n+1:2*n)=0.0 ! bulk
phi(0) = phi(1) ! reflection at x = 0

! Recover xtotal from input
xtotal(:,n+1:2*n)=0.0 ! bulk
do is = 1,Npoorsv
   xtotal(is,1:n) = x(1+n*is:n+n*is) !  segments volume fractions (xtotal) are read from the x provided by kinsol
enddo


!inverse_of_vpolvsol=1/(vpol*vsol)

avpos=0.0
avneg=0.0

do j=1,ntot
   avpos(j)=vpos*expmupos*xh(j)**vpos*dexp(-phi(j)) ! volume fraction of cations
   avneg(j)=vneg*expmuneg*xh(j)**vneg*dexp(phi(j)) ! volume fraction of anions
   avHplus(j)=expmuHplus*xh(j)*dexp(-phi(j)) ! volume fraction of H+
   avOHmin(j)=expmuOHmin*xh(j)*dexp(phi(j)) ! volume fraction of OH-

   if (Nacids.gt.0) then
     do is=1,Nacids
       fAmin(is,j)=1.0/(avHplus(j)/(Ka(is)*xh(j))+1.0)
     enddo
   endif

   if (Nbasics.gt.0) then
     do is=1,Nbasics
       fBHplus(is,j)=1.0/(avOHmin(j)/(Kb(is)*xh(j))+1.0)
     enddo
   endif

enddo



! Caculation of dielectric function
! Everything that it is not water or ions has dielectric dielP

do j = 1, ntot
dielpol(j) = 1.0 - xh(j) - avpos(j) - avneg(j) - avHplus(j) - avOHmin(j)
enddo 
call dielectfcn(dielpol,epsfcn,Depsfcn)



! Calculation of xpot
do i = 1, ntot
! osmotic pressure
   xpot(0,i) = xh(i)**vpol(0) ! exp(-pi(r)v_pol) / units of v_pol: nm^3
! dielectrics
   gradphi2 = ((phi(i+1)-phi(i))/delta)**2
   xpot(0,i) = xpot(0,i)*exp(Depsfcn(i)*gradphi2*vpol(0)*vsol*wperm/2.0)
enddo 



do i = 1, ntot
  do is = 1, Npoorsv

!   calculate xpot(i, is)

    protemp = 0.0
    
      do js = 1, Npoorsv 
         do j = 1, ntot
            protemp = protemp+st(is,js)*Xu(i,j,is,js)*xtotal(js,j)/(vpol(js)*vsol) ! vpol*vsol in units of nm^3
!           protemp = protemp+st(is,js)/(vpol*vsol)*Xu(i,j,is,js)*xtotal(js,j) ! vpol*vsol in units of nm^3
         enddo
      enddo

    xpot(is,i) = xh(i)**vpol(is) ! exp(-pi(r)v_pol) / units of v_pol: nm^3
! dielectrics
    gradphi2 = ((phi(i+1)-phi(i))/delta)**2
    xpot(is,i) = xpot(is,i)*exp(Depsfcn(i)*gradphi2*vpol(is)*vsol*wperm/2.0)
    xpot(is,i) = xpot(is,i)*dexp(protemp) 

  enddo

!  do ic = 1,Ncharge
!    protemp=-phi(i)*float(charge(ic))
!    xpotc(ic,i)=dexp(protemp)
!  enddo

  if (Nacids.gt.0) then
    do ic = 1,Nacids
      xpot_a(ic,i) = 1.0/fAmin(ic,i)*exp(phi(i))
    enddo
  endif

  if (Nbasics.gt.0) then
    do ic = 1,Nbasics
      xpot_b(ic,i)= 1.0/fBHplus(ic,i)*exp(-phi(i))
    enddo
  endif

enddo


do is = 0,Npoorsv
   xpot(is,n+1:2*n)=xpot(is,n)
enddo

if (Nacids.gt.0) then
  do ic = 1,Nacids
    xpot_a(ic,n+1:2*n)=xpot_a(ic,n)
  enddo
endif

if (Nbasics.gt.0) then
  do ic = 1,Nbasics
    xpot_b(ic,n+1:2*n)=xpot_b(ic,n)
  enddo
endif

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


do ii=1,maxntotcounter ! position of center of mass 

   do i=1, cuantas ! loop over conformations
 
      pro(i) = exp(-Uchain(i))

      do j=minpos(i,ii), maxpos(i,ii) ! loop over lattice position

         k = j-minpos(i,ii)+1 ! k may be lager than ntot

         do is = 0, Npoorsv 
            pro(i)= pro(i) * xpot(is,j)**inn(is,i,ii,k)
         enddo

         if (Nacids.gt.0) then
           do ic = 1, Nacids
             pro(i)= pro(i) * xpot_a(ic,j)**inn_a(ic,i,ii,k) 
           enddo
         endif
     
         if (Nbasics.gt.0) then
           do ic = 1, Nbasics
             pro(i)= pro(i) * xpot_b(ic,j)**inn_b(ic,i,ii,k)
           enddo
         endif
 
      enddo !j

      all_tosend(ii) = all_tosend(ii) + pro(i) !q_tosend(ii)=q_tosend(ii)+pro(i)
      all_tosend(ntot+ii) = all_tosend(ntot+ii) + pro(i)*dlog(pro(i)) !sumprolnpro_tosend(ii) = sumprolnpro_tosend(ii) + pro(i)*dlog(pro(i))
      all_tosend(ntot*2+ii) = all_tosend(ntot*2+ii) + pro(i)*Uchain(i) ! sumprouchain_tosend(ii) = sumprouchain_tosend(ii) + pro(i)*Uchain(i)
      all_tosend(ntot*3+ii) = all_tosend(ntot*3+ii) + pro(i) ! xpol_tosend(ii) = xpol_tosend(ii)+pro(i)


      do j = 1, long ! loop over number of segments
         sumtrans_tosend(ii,j) =  sumtrans_tosend(ii,j) +  pro(i)*float(Ntrans(j,i))
      enddo

      do j=minpos(i,ii), maxpos(i,ii)
         k = j-minpos(i,ii)+1 ! k may be larger than ntot

         do is = 0, Npoorsv 
            avpol_tmp(is,j)=avpol_tmp(is,j)+pro(i)*inn(is,i,ii,k)*factorcurv(ii,j) ! avpol_tmp is avg number of segments "is" at position "j" 
         enddo

         if (Nacids.gt.0) then
           do ic = 1,Nacids
             avpola_tmp(ic,j)=avpola_tmp(ic,j)+pro(i)*inn_a(ic,i,ii,k)*factorcurv(ii,j) ! avpola_tmp is avg number of acid segments "ic" at position "j"
           enddo
         endif

         if (Nbasics.gt.0) then
           do ic = 1,Nbasics
             avpolb_tmp(ic,j)=avpolb_tmp(ic,j)+pro(i)*inn_b(ic,i,ii,k)*factorcurv(ii,j) ! avpolb_tmp is avg number of basic segments "ic" at position "j" 
           enddo
         endif

      enddo

   enddo ! i

enddo   ! ii

avpol_tosend(:, 1:ntot)=avpol_tmp(:, 1:ntot) 
avpola_tosend(:, 1:ntot)=avpola_tmp(:, 1:ntot)
avpolb_tosend(:, 1:ntot)=avpolb_tmp(:, 1:ntot)

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
!   call MPI_ALLREDUCE(xpol_tosend, xpol, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(all_tosend, all_toreceive, ntot*4, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD, err)
!   call MPI_ALLREDUCE(q_tosend, q, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
!   call MPI_ALLREDUCE(sumprolnpro_tosend, sumprolnpro, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
!   call MPI_ALLREDUCE(sumprouchain_tosend, sumprouchain, ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
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

q(:) = all_toreceive(1:ntot)
sumprolnpro(:) = all_toreceive(ntot+1:2*ntot)
sumprouchain(:) = all_toreceive(ntot*2+1:3*ntot)
xpol(:) = all_toreceive(ntot*3+1:4*ntot)


sumpol = 0.0

do i = 1, ntot
   do is = 0, Npoorsv
 
      select case (curvature)

       case (0)
        sumpol = sumpol + avpol(is,i)*delta ! final result in units of chains/nm^2
       case(1)
        sumpol = sumpol + avpol(is,i)*(float(i)-0.5)*delta*delta*2.0*pi ! final result in units of chains/nm
       case(2)
        sumpol = sumpol + avpol(is,i)*(((float(i)-0.5)*delta)**2)*delta*4.0*pi ! final result in units of chains/micelle

      end select

   enddo
enddo


sumpol = sumpol/(vchain*vsol) 
avpol = avpol/sumpol*npol ! integral of avpol is fixed
avpola = avpola/sumpol*npol
avpolb = avpolb/sumpol*npol
sumpol = 0.0


do i = 1, ntot
   select case (curvature)
    case (0)
     sumpol = sumpol + xpol(i)*delta ! final result in units of chains/nm^2
    case(1)
     sumpol = sumpol + xpol(i)*(float(i)-0.5)*delta*delta*2.0*pi ! final result in units of chains/nm
    case(2)
     sumpol = sumpol + xpol(i)*(((float(i)-0.5)*delta)**2)*delta*4.0*pi ! final result in units of chains/micelle
   end select
enddo

xpol = xpol/sumpol*npol ! integral of avpol is fixed

trans = 0.0

do i = 1, maxntotcounter
   select case (curvature)
    case (0)
     trans(:) = trans(:) + sumtrans(i,:)/q(i)*xpol(i)*delta ! final result in units of chains/nm^2
    case(1)
     trans(:) = trans(:) + sumtrans(i,:)/q(i)*xpol(i)*(float(i)-0.5)*delta*delta*2.0*pi ! final result in units of chains/nm
    case(2)
     trans(:) = trans(:) + sumtrans(i,:)/q(i)*xpol(i)*(((float(i)-0.5)*delta)**2)*delta*4.0*pi ! final result in units of chains/micelle
   end select
enddo

trans(:) = trans(:)/npol

! contruction of f and the volume fractions

do i=1,n

   f(i)=xh(i)+avneg(i)+avpos(i)+avHplus(i)+avOHmin(i)-1.0d0

   do is=0, Npoorsv
      f(i) = f(i) + avpol(is,i)
   enddo

enddo

xcharge(:) = avpos(:)/(vpos*vsol)-avneg(:)/(vneg*vsol)+avHplus(:)/vsol-avOHmin(:)/vsol ! xcharge is avg charge density

inverse_of_two=1/2.0

do i = 1,ntot

   if (Nacids.gt.0) then
     do ic= 1,Nacids
       xcharge(i)=xcharge(i)-avpola(ic,i)*fAmin(ic,i)/(vpol_a(ic)*vsol)
     enddo
   endif

   if (Nbasics.gt.0) then
     do ic= 1,Nbasics
       xcharge(i)=xcharge(i)+avpolb(ic,i)*fBHplus(ic,i)/(vpol_b(ic)*vsol)
     enddo
   endif

   select case (curvature)
    case (0)

!     f(i+n*(Npoorsv+1))=xcharge(i)*inverse_of_wperm+(phi(i+1)-2*phi(i)+phi(i-1))*delta**(-2)  


     f(i+n*(Npoorsv+1))=xcharge(i) + wperm*epsfcn(i)*(phi(i+1)-2*phi(i)+phi(i-1))*delta**(-2) &
     +wperm*(epsfcn(i+1)-epsfcn(i))*(phi(i+1)-phi(i))*delta**(-2)

    case(1)
     f(i+n*(Npoorsv+1))=xcharge(i) + wperm*epsfcn(i)*(phi(i+1)-phi(i))*delta**(-2)/(float(i)-0.5) &
     +wperm*epsfcn(i)*(phi(i+1)-2.0*phi(i)+phi(i-1))*delta**(-2) &
     +wperm*(epsfcn(i+1)-epsfcn(i))*(phi(i+1)-phi(i))*delta**(-2)  

    case(2)
     f(i+n*(Npoorsv+1))=xcharge(i) + 2.0*wperm*epsfcn(i)*(phi(i+1)-phi(i))*delta**(-2)/(float(i)-0.5) &
     +wperm*epsfcn(i)*(phi(i+1)-2.0*phi(i)+phi(i-1))*delta**(-2) &
     +wperm*(epsfcn(i+1)-epsfcn(i))*(phi(i+1)-phi(i))*delta**(-2)

    end select

   f(i+n*(Npoorsv+1))=-f(i+n*(Npoorsv+1))*inverse_of_two
!   f(i+n*(Npoorsv+1))=f(i+n*(Npoorsv+1))/(-2.0)

enddo

do is=1,Npoorsv
   do i=1,n ! xtotal
      f(i+n*is) = -avpol(is,i)+xtotal(is,i) 
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
   print*, iter, algo, q(1), Q(2), q(3), q(4)
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

