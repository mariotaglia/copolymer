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
integer NC,lastlayer
real*8 all_tosend(4*ntot), all_toreceive(4*ntot)
integer*4 ier2
real*8 protemp
real*8 x((Npoorsv+2)*ntot),f((Npoorsv+2)*ntot)
real*8 xh(dimR+1,dimZ) 
real*8 xpot(0:Npoorsv,dimR+20,dimZ), xpot_a(0:Nacids,2*dimR,dimZ), xpot_b(0:Nbasics,2*dimR,dimZ)
real*8,allocatable :: pro(:)
!real*8 time1, time2, duration, looptime1, looptime2, loopduration
integer iR,iZ,kZ,kkZ,k,i,j,ic,aR,aZ,iZm,iZp,jZp,jZm        ! dummy indices
integer is, js,ia,ib,iiR,iiZ,jR,jZ
integer err
integer n
real*8 avpol_tmp(0:Npoorsv,2*dimR,dimZ), avpola_tmp(0:Nacids,2*dimR,dimZ), avpolb_tmp(0:Nbasics,2*dimR,dimZ) ! overdim R coordinate just in case
real*8 avpol_tosend(0:Npoorsv,dimR,dimZ), avpola_tosend(0:Nacids,dimR,dimZ), avpolb_tosend(0:Nbasics,dimR,dimZ)
real*8 xpol_tosend(dimR,dimZ)
real*8 algo, algo1,algo2
double precision, external :: factorcurv
real*8 sumpol
real*8 q_tosend(dimR,dimZ)
real*8 sumprolnpro_tosend(dimR,dimZ), sumprouchain_tosend(dimR,dimZ)
real*8 sumtrans_tosend(dimR,dimZ,maxlong)
real*8 sumtrans(dimR,dimZ,maxlong)
real*8 gradphi2
integer, external :: PBCSYMI
integer, external :: PBCREFI
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

xh(dimR+1,:)=1.0


do iR=1,dimR
do iZ=1,dimZ
  phi(iR,iZ)=x(n*(Npoorsv+1)+dimR*(iZ-1)+iR)
enddo
enddo

phi(0,:)=phi(1,:) ! symmetry at r = 0
phi(dimR+1,:)=0.0 ! bulk for r -> inf

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
   jZp=iZ+1 ! jZ plus one
   jZm=iZ-1 ! jZ minus one

   if(PBCflag.eq.1) then
      iZp=PBCSYMI(jZp,dimZ)
      iZm=PBCSYMI(jZm,dimZ)
   else if(PBCflag.eq.2) then
      iZp=PBCREFI(jZp,dimZ)
      iZm=PBCREFI(jZm,dimZ)
   endif

do iR = 1,dimR
   gradphi2 = ((phi(iR+1,iZ)-phi(iR-1,iZ))/2.0/deltaR)**2+((phi(iR,iZp)-phi(iR,iZm))/2.0/deltaZ)**2
   xpot(0,iR,iZ) = xpot(0,iR,iZ)*exp(Depsfcn(iR,iZ)*gradphi2*vpol(0)*vsol*wperm/2.0)
enddo 
enddo

do iZ = 1, dimZ
do iR = 1, dimR

  do is = 1, Npoorsv
!   calculate xpot(i, is)

    protemp = 0.0
    
      do js = 1, Npoorsv 
         do jR = Rini_kais, Rfin_kais
         do jZ = -Xulimit, Xulimit
            kZ=jZ+iZ

            if(PBCflag.eq.1)kkZ=PBCSYMI(kZ,dimZ)
            if(PBCflag.eq.2)kkZ=PBCREFI(kZ,dimZ)

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
  jZm=iZ-1

   if(PBCflag.eq.1) then
      iZp=PBCSYMI(jZp,dimZ)
      iZm=PBCSYMI(jZm,dimZ)
   else if(PBCflag.eq.2) then
      iZp=PBCREFI(jZp,dimZ)
      iZm=PBCREFI(jZm,dimZ)
   endif


    do iR= 1, dimR
      gradphi2 = ((phi(iR+1,iZ)-phi(iR-1,iZ))/2.0/deltaR)**2+((phi(iR,iZp)-phi(iR,iZm))/2.0/deltaZ)**2
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    probability distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

avpola = 0.0
avpolb = 0.0
avpol = 0.0
xpol = 0.0
q = 0.0
sumprolnpro = 0.0
sumprouchain=0.0
sumtrans_tosend = 0.0
sumtrans = 0.0
all_toreceive = 0.0
xpot_a(0,:,:)=1.0
xpot_b(0,:,:)=1.0

do NC = 1, Ncomp ! loop over components

allocate(pro(cuantas(NC)))

q_tosend=0.0d0                   ! init q to zero for each component
avpola_tosend = 0.0
avpolb_tosend = 0.0
avpol_tosend = 0.0
xpol_tosend = 0.0
sumprolnpro_tosend = 0.0
sumprouchain_tosend=0.0
all_tosend = 0.0
avpola_tmp = 0.0
avpolb_tmp = 0.0
avpol_tmp = 0.0

do iiR=minntotR(NC), maxntotR(NC) ! position of center of mass 
do iiZ=minntotZ(NC), maxntotZ(NC)
 
   do i=1, cuantas(NC) ! loop over conformations

      pro(i) = exp(-Uchain(i, NC))

      do k=1,long(NC)
         iZ = innZ(k,i,NC)+iiZ

         if(PBCflag.eq.1)aZ = PBCSYMI(iZ,dimZ)
         if(PBCflag.eq.2)aZ = PBCREFI(iZ,dimZ)

         aR = innR(k,i,iiR,NC)
         is = segpoorsv(k,NC)
         ia = acidtype(k,NC)
         ib = basictype(k,NC) 
                  
         if (aR.gt.dimR) then
           pro(i) = pro(i) * xpot_bulk(is)
           pro(i) = pro(i) * xpota_bulk(ia)
           pro(i) = pro(i) * xpotb_bulk(ib)
         else
           pro(i)= pro(i) * xpot(is,aR,aZ)
           pro(i)= pro(i) * xpot_a(ia,aR,aZ) 
           pro(i)= pro(i) * xpot_b(ib,aR,aZ)
         endif

      enddo !k

      q_tosend(iiR,iiZ) = q_tosend(iiR,iiZ) + pro(i) ! all_tosend(ii) = all_tosend(ii) + pro(i) 
      sumprolnpro_tosend(iiR,iiZ) = sumprolnpro_tosend(iiR,iiZ) + pro(i)*dlog(pro(i)) ! all_tosend(ntot+ii) = all_tosend(ntot+ii) + pro(i)*dlog(pro(i))
      sumprouchain_tosend(iiR,iiZ) = sumprouchain_tosend(iiR,iiZ) + pro(i)*Uchain(i,NC) !all_tosend(ntot*2+ii) = all_tosend(ntot*2+ii) + pro(i)*Uchain(i) 
      xpol_tosend(iiR,iiZ) = xpol_tosend(iiR,iiZ)+pro(i) ! all_tosend(ntot*3+ii) = all_tosend(ntot*3+ii) + pro(i)
      
      do j = 1, long(NC) ! loop over number of segments

         iZ = innZ(j,i,NC)+iiZ
         if(PBCflag.eq.1)aZ = PBCSYMI(iZ,dimZ)
         if(PBCflag.eq.2)aZ = PBCREFI(iZ,dimZ)
         
         aR = innR(j,i,iiR,NC)
         is = segpoorsv (j,NC)
         ia = acidtype (j,NC)
         ib = basictype (j,NC)

         sumtrans_tosend(iiR,iiZ,j) =  sumtrans_tosend(iiR,iiZ,j) +  pro(i)*float(Ntrans(j,i,NC))

         avpol_tmp(is,aR,aZ) = avpol_tmp(is,aR,aZ)+pro(i)*factorcurv(iiR,aR) ! avpol_tmp is avg number of segments "is" at position "j" 
         avpola_tmp(ia,aR,aZ) = avpola_tmp(ia,aR,aZ)+pro(i)*factorcurv(iiR,aR) ! avpola_tmp is avg number of acid segments "ic" at position "j"
         avpolb_tmp(ib,aR,aZ) = avpolb_tmp(ib,aR,aZ)+pro(i)*factorcurv(iiR,aR) ! avpolb_tmp is avg number of basic segments "ic" at position "j" 
      enddo ! j

   enddo ! i

enddo ! iiR
enddo ! iiZ

!! sum of sements from GC chains in the bulk

if (flagGC(NC).eq.1) then

  do iiR = maxntotR(NC)+1,maxntotR(NC)+30
  do iiZ = minntotZ(NC),maxntotZ(NC)
    do i = 1,cuantas(NC)
      do j = 1, long(NC) ! loop over number of segments

         iZ = innZ(j,i,NC)+iiZ
         if(PBCflag.eq.1)aZ = PBCSYMI(iZ,dimZ)
         if(PBCflag.eq.2)aZ = PBCREFI(iZ,dimZ)

         aR = innR(j,i,iiR,NC)
         is = segpoorsv (j,NC)
         ia = acidtype (j,NC)
         ib = basictype (j,NC)

         avpol_tmp(is,aR,aZ) = avpol_tmp(is,aR,aZ) + probulk(NC) 
         avpola_tmp(ia,aR,aZ) = avpola_tmp(ia,aR,aZ) + probulk(NC) 
         avpolb_tmp(ib,aR,aZ) = avpolb_tmp(ib,aR,aZ) + probulk(NC)

      enddo !j
    enddo !i 
  enddo !iiZ
  enddo !iiR

endif 


avpol_tosend(:, 1:dimR, 1:dimZ)=avpol_tmp(:, 1:dimR, 1:dimZ) 
avpola_tosend(:, 1:dimR, 1:dimZ)=avpola_tmp(:, 1:dimR, 1:dimZ)
avpolb_tosend(:, 1:dimR, 1:dimZ)=avpolb_tmp(:, 1:dimR, 1:dimZ)


!------------------ MPI -----------------`-----------------------------

!call MPI_Barrier(MPI_COMM_WORLD, err)

   call MPI_ALLREDUCE(avpola_tosend, avpola(:,:,:,NC), (Nacids+1)*ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(avpolb_tosend, avpolb(:,:,:,NC), (Nbasics+1)*ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(avpol_tosend, avpol(:,:,:,NC), (Npoorsv+1)*ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(xpol_tosend, xpol(:,:,NC), ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(q_tosend, q(:,:,NC), ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(sumprolnpro_tosend, sumprolnpro(:,:,NC), ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(sumprouchain_tosend, sumprouchain(:,:,NC), ntot, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)
   call MPI_ALLREDUCE(sumtrans_tosend, sumtrans, ntot*maxlong, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, err)

!----------------------- Norm -----------------------------------------


sumpol = 0.0

do iR = 1, dimR
do iZ = 1, dimZ
do is = 0, Npoorsv
 select case (curvature)
  case (0)
   sumpol = sumpol + avpol(is,iR,iZ,NC)*deltaR*deltaZ ! final result in units of chains/nm^2 (1D) or in units of chains/nm of belt (2D)
  case(1)
   sumpol = sumpol + avpol(is,iR,iZ,NC)*deltaZ*(float(iR+dimRini)-0.5)*deltaR*deltaR*2.0*pi ! final result in units of chains/nm of fiber (1D) or chains/fiber (2D)
  case(2)
   sumpol = sumpol + avpol(is,iR,iZ,NC)*(((float(iR+dimRini)-0.5)*deltaR)**2)*deltaR*4.0*pi ! final result in units of chains/micelle
 end select
enddo
enddo
enddo

if (flagGC(NC).eq.0) then
   sumpol = sumpol/(vchain(NC)*vsol) 
   avpol(:,:,:,NC) = avpol(:,:,:,NC)/sumpol*npol*npolratio(NC) ! integral of avpol is fixed
   avpola(:,:,:,NC) = avpola(:,:,:,NC)/sumpol*npol*npolratio(NC)
   avpolb(:,:,:,NC) = avpolb(:,:,:,NC)/sumpol*npol*npolratio(NC)
else
   do is=0,Npoorsv
      avpol(is,:,:,NC) = avpol(is,:,:,NC)*expmupol(NC)*vpol(is) ! avpol = rho_pol * vpol = q * expmupol * vpol / vsol
   enddo
   do ia=1,Nacids
      avpola(ia,:,:,NC) = avpola(ia,:,:,NC)*expmupol(NC)*vpol_a(ia)
   enddo
   do ib=1,Nbasics
      avpolb(ib,:,:,NC) = avpolb(ib,:,:,NC)*expmupol(NC)*vpol_b(ib)
   enddo
endif

sumpol = 0.0

do iR = 1, dimR
do iZ = 1, dimZ
  select case (curvature)
   case (0)
    sumpol = sumpol + xpol(iR,iZ,NC)*deltaR*deltaZ ! final result in units of chains/nm^2 (1D) or chains/nm of belt (2D)
   case (1)
    sumpol = sumpol + xpol(iR,iZ,NC)*deltaZ*(float(iR+dimRini)-0.5)*deltaR*deltaR*2.0*pi ! final result in units of chains/nm (1D) or chains/fiber (2D)
   case (2)
    sumpol = sumpol + xpol(iR,iZ,NC)*(((float(iR+dimRini)-0.5)*deltaR)**2)*deltaR*4.0*pi ! final result in units of chains/micelle
   end select
enddo
enddo

if (flagGC(NC).eq.0) xpol(:,:,NC) = xpol(:,:,NC)/sumpol*npol*npolratio(NC) ! integral of avpol is fixed 
if (flagGC(NC).eq.1) xpol(:,:,NC) = xpol(:,:,NC)*expmupol(NC)/vsol

trans(:,NC) = 0.0

do iR = minntotR(NC), maxntotR(NC)
do iZ = minntotZ(NC), maxntotZ(NC)
  select case (curvature)
   case (0)
    trans(:,NC) = trans(:,NC) + sumtrans(iR,iZ,:)/q(iR,iZ,NC)*xpol(iR,iZ,NC)*deltaR*deltaZ ! final result in units of chains/nm^2 (1D) or chains/nm of belt (2D)
   case(1)
    trans(:,NC) = trans(:,NC) + sumtrans(iR,iZ,:)/q(iR,iZ,NC)*xpol(iR,iZ,NC)*deltaZ*(float(iR+dimRini)-0.5)*deltaR*deltaR*2.0*pi ! final result in units of chains/nm (1D) or chains/fiber (2D)
   case(2)
    trans(:,NC) = trans(:,NC) + sumtrans(iR,iZ,:)/q(iR,iZ,NC)*xpol(iR,iZ,NC)*(((float(iR+dimRini)-0.5)*deltaR)**2)*deltaR*4.0*pi ! final result in units of chains/micelle
   end select
enddo
enddo

if (flagGC(NC).eq.0) trans(:,NC) = trans(:,NC)/npol/npolratio(NC)
if (flagGC(NC).eq.1) trans(:,NC) = trans(:,NC)*expmupol(NC)/vsol

deallocate(pro)

! if (flagGC(NC).eq.1) print*,sumpol*expmupol(NC)/totalcuantas(NC), expmupol(NC)/totalcuantas(NC)

      !!!!!!
enddo ! NC !
      !!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contruction of f and the volume fractions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do iR = 1, dimR
do iZ = 1, dimZ

   f(dimR*(iZ-1)+iR)=xh(iR,iZ)+avneg(iR,iZ)+avpos(iR,iZ)+avHplus(iR,iZ)+avOHmin(iR,iZ)-1.0d0

   do NC = 1,Ncomp
   do is=0, Npoorsv
      f(dimR*(iZ-1)+iR) = f(dimR*(iZ-1)+iR) + avpol(is,iR,iZ,NC)
   enddo
   enddo

enddo
enddo

xcharge(:,:) = avpos(:,:)/(vpos*vsol)-avneg(:,:)/(vneg*vsol)+avHplus(:,:)/vsol-avOHmin(:,:)/vsol ! xcharge is avg charge density


do iR = 1, dimR
do iZ = 1, dimZ

   do NC = 1, Ncomp
   do ic= 1,Nacids
     xcharge(iR,iZ)=xcharge(iR,iZ)-avpola(ic,iR,iZ,NC)*fAmin(ic,iR,iZ)/(vpol_a(ic)*vsol)
   enddo

   do ic= 1,Nbasics
     xcharge(iR,iZ)=xcharge(iR,iZ)+avpolb(ic,iR,iZ,NC)*fBHplus(ic,iR,iZ)/(vpol_b(ic)*vsol)
   enddo
   enddo ! NC

   iZp=iZ+1
   iZm=iZ-1

   if(PBCflag.eq.1) then
   jZp = PBCSYMI(iZp,dimZ)
   jZm = PBCSYMI(iZm,dimZ)
   else if (PBCflag.eq.2) then
   jZp = PBCREFI(iZp,dimZ)
   jZm = PBCREFI(iZm,dimZ)
   endif 

   select case (curvature)

    case (0)
     f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=xcharge(iR,iZ) &
     +wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-2.0*phi(iR,iZ)+phi(iR-1,iZ))*deltaR**(-2) &
     +wperm*epsfcn(iR,iZ)*(phi(iR,jZp)-2.0*phi(iR,iZ)+phi(iR,jZm))*deltaZ**(-2) &
     +wperm*(epsfcn(iR+1,iZ)-epsfcn(iR-1,iZ))/2.0*(phi(iR+1,iZ)-phi(iR-1,iZ))/2.0*deltaR**(-2) &
     +wperm*(epsfcn(iR,jZp)-epsfcn(iR,jZm))*(phi(iR,jZp)-phi(iR,jZm))/4.0*deltaZ**(-2) 

    case(1)
     f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=xcharge(iR,iZ) &
     + wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2)/(float(iR+dimRini)-0.5) &
     +wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-2.0*phi(iR,iZ)+phi(iR-1,iZ))*deltaR**(-2) &
     +wperm*epsfcn(iR,iZ)*(phi(iR,jZp)-2.0*phi(iR,iZ)+phi(iR,jZm))*deltaZ**(-2) &
     +wperm*(epsfcn(iR+1,iZ)-epsfcn(iR,iZ))*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2) &
     +wperm*(epsfcn(iR,jZp)-epsfcn(iR,jZm))*(phi(iR,jZp)-phi(iR,jZm))/4.0*deltaZ**(-2)

    case(2)
     f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=xcharge(iR,iZ) &
     + 2.0*wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2)/(float(iR+dimRini)-0.5) &
     +wperm*epsfcn(iR,iZ)*(phi(iR+1,iZ)-2.0*phi(iR,iZ)+phi(iR-1,iZ))*deltaR**(-2) &
     +wperm*(epsfcn(iR+1,iZ)-epsfcn(iR,iZ))*(phi(iR+1,iZ)-phi(iR,iZ))*deltaR**(-2)

    end select
   
   f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)=f(n*(Npoorsv+1)+dimR*(iZ-1)+iR)/(-2.0)

enddo
enddo

do is=1,Npoorsv
   do iR=1,dimR
   do iZ=1,dimZ ! xtotal
      f(n*is+dimR*(iZ-1)+iR) = xtotal(is,iR,iZ) 
      do NC = 1,Ncomp
          f(n*is+dimR*(iZ-1)+iR) = -avpol(is,iR,iZ,NC)+f(n*is+dimR*(iZ-1)+iR) 
      enddo ! NC
   enddo
   enddo
enddo ! is


iter=iter+1

algo = 0.0
algo1 = 0.0
algo2 = 0.0

do i = 1, n*(Npoorsv+2)
   algo = algo + f(i)**2
end do

if(rank.eq.0) then 
   print*, iter, algo, q(1,1,1)
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

