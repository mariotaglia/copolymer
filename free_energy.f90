subroutine calc_free_energy

use globals
use mcharge
use partfunc
use pis
use MPI
use bulk
use layer
use volume
use mkai
implicit none                                 

real*8 gradphi2
real*8 Free_energy, F_Mix_s, F_Mix_pos, F_mix_p
real*8 F_parCopA, F_parMol, F_Asoc
real*8 F_Mix_neg, F_Mix_Hplus                 
real*8 Free_energy2, sumpi, sumrho, sumel,sumpol, mupol, sumdiel
real*8 sumeq ! LEO
real*8 F_Mix_OHmin, F_Conf, F_Uchain     
real*8  F_Eq,F_vdW(Npoorsv,Npoorsv), F_electro                            
integer is, js, iR, iZ, jR, jZ, kZ, kkZ, jZp, iZp, iZm, jZm
integer, external :: PBCSYMI
integer, external :: PBCREFI
double precision, external :: jacobian
integer NC,MC
real*8 cteAc, cteBasic !LEO
character*17 F_vdWfilename(Npoorsv,Npoorsv)

if(rank.eq.0) then 
print*, 'Starting free energy calculation...'
! open files
open(unit=301, file='F_tot.dat')                           
open(unit=302, file='F_mixs.dat')                          
open(unit=303, file='F_mixpos.dat')                        
open(unit=304, file='F_mixneg.dat')                        
open(unit=305, file='F_mixH.dat')                          
open(unit=306, file='F_mixOH.dat')                         
open(unit=307, file='F_conf.dat')                          
open(unit=308, file='F_eq.dat')                            

do is=1,Npoorsv
do js=1,Npoorsv
write(F_vdWfilename(is,js),'(A6,BZ,I3.3,A1,I3.3,A4)')'F_vdW.',is,'.',js,'.dat'
open(unit=10000*is+js, file=F_vdWfilename(is,js) )                           
enddo
enddo
open(unit=3110, file='F_electro.dat')                       
open(unit=312, file='F_tot2.dat')                          
open(unit=313, file='F_mixp.dat')                          
open(unit=314, file='F_Uchain.dat')
endif


!---------------------------------------------------------------------
!        Calculate Free Energy                                        
!-------------------------------------------------------------------- 

Free_Energy = 0.0                                                
Free_Energy2 = 0.0                                               


! 1. Solvent translational entropy
F_Mix_s = 0.0                                                    

do iR = 1, dimR
do iZ = 1, dimZ
   F_Mix_s = F_Mix_s + xsol(iR,iZ)*(dlog(xsol(iR,iZ))-1.0)*jacobian(iR)*deltaR*deltaZ/vsol ! mix entropy of solvent                                      
   F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)*jacobian(iR)*deltaR*deltaZ/vsol ! with respect to bulk                             
enddo                                                             
enddo

Free_Energy = Free_Energy + F_Mix_s                              

!F_mup = 0.0                                                    
!do iC = 1, maxntotcounter                                                
!F_mup = F_mup - xpol(iC)*jacobian(iC)*delta*mupol                                      
!enddo                                                            
!Free_Energy = Free_Energy + F_mup                              

! 1. Polymer translational entropy
F_Mix_p = 0.0                                                    

do NC = 1, Ncomp
do iR = 1, maxntotcounterR                                                
do iZ = 1, maxntotcounterZ
   F_Mix_p = F_Mix_p + xpol(iR,iZ,NC)*(dlog(xpol(iR,iZ,NC))-1.0)*jacobian(iR)*deltaR*deltaZ ! mix entropy of chains with respect to bulk (xpolbulk=0) 
enddo                                                            
enddo
enddo

Free_Energy = Free_Energy + F_Mix_p            

! 2. Pos ion translational entropy

F_Mix_pos = 0.0                                                  

do iR=1,dimR
do iZ=1,dimZ
  F_Mix_pos = F_Mix_pos + avpos(iR,iZ)*(dlog(avpos(iR,iZ)/vpos)-1.0-dlog(expmupos))*jacobian(iR)*deltaR*deltaZ/(vsol*vpos) ! mix entropy of cations and N*Mu term
  F_Mix_pos = F_Mix_pos - xposbulk*(dlog(xposbulk/vpos)-1.0-dlog(expmupos))*jacobian(iR)*deltaR*deltaZ/(vsol*vpos) ! with respect to bulk
enddo
enddo 

Free_energy = Free_Energy + F_Mix_pos

!3. Neg ion translational entropy

F_Mix_neg = 0.0                                                  

do iR = 1, dimR
do iZ = 1, dimZ
  F_Mix_neg = F_Mix_neg + avneg(iR,iZ)*(dlog(avneg(iR,iZ)/vpos)-1.0-dlog(expmuneg))*jacobian(iR)*deltaR*deltaZ/(vsol*vneg) ! mix entropy and N*Mu term of anions 
  F_Mix_neg = F_Mix_neg - xnegbulk*(dlog(xnegbulk/vpos)-1.0-dlog(expmuneg))*jacobian(iR)*deltaR*deltaZ/(vsol*vneg) ! with respect to bulk
enddo
enddo

Free_Energy = Free_Energy + F_Mix_neg                            

! 4. H+ translational entropy

F_Mix_Hplus = 0.0                                                

do iR = 1, dimR
do iZ = 1, dimZ
  F_Mix_Hplus = F_Mix_Hplus + avHplus(iR,iZ)/vsol*(dlog(avHplus(iR,iZ))-1.0-dlog(expmuHplus))*jacobian(iR)*deltaR*deltaZ ! mix entropy and N*Mu term of Hplus
  F_Mix_Hplus = F_Mix_Hplus - xHplusbulk/vsol*(dlog(xHplusbulk)-1.0-dlog(expmuHplus))*jacobian(iR)*deltaR*deltaZ ! with respect to bulk      
enddo                                                            
enddo

Free_Energy = Free_Energy + F_Mix_Hplus                          

! 5. OH- translational entropy 

F_Mix_OHmin = 0.0                                                

do iR = 1, dimR
do iZ = 1, dimZ
  F_Mix_OHmin = F_Mix_OHmin + avOHmin(iR,iZ)/vsol*(dlog(avOHmin(iR,iZ))-1.0-dlog(expmuOHmin))*jacobian(iR)*deltaR*deltaZ ! mix entropy and N*Mu term of OHmin          
  F_Mix_OHmin = F_Mix_OHmin - xOHminbulk/vsol*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*jacobian(iR)*deltaR*deltaZ ! with respect to bulk       
enddo                                                            
enddo

Free_Energy = Free_Energy + F_Mix_OHmin                          

! 6. Polymer conformational entropy                                         

do NC = 1, Ncomp
F_conf = 0 

do iR = 1, maxntotcounterR
do iZ = 1, maxntotcounterZ
  F_Conf = F_conf + (sumprolnpro(iR,iZ,NC)/q(iR,iZ,NC)-dlog(q(iR,iZ,NC)))*jacobian(iR)*deltaR*deltaZ*xpol(iR,iZ,NC)
enddo 
enddo

Free_Energy = Free_Energy + F_Conf

F_Uchain = 0.0

do iR=1, maxntotcounterR
do iZ=1, maxntotcounterZ
  F_Uchain = F_Uchain + deltaR*deltaZ*xpol(iR,iZ,NC)*jacobian(iR)*(sumprouchain(iR,iZ,NC)/q(iR,iZ,NC))
enddo
enddo

Free_Energy = Free_Energy + F_Uchain
enddo

! 7. Chemical Equilibrium                                              

F_Eq = 0.0                                                       

!do NC = 1, Ncomp

do iR=1,dimR
cteAc = jacobian(iR)*deltaR*deltaZ/(vpol_a(1)*vsol)
cteBasic = jacobian(iR)*deltaR*deltaZ/(vpol_b(1)*vsol)
do iZ=1,dimZ

!LEO

!eq acido base del acido:
F_Eq = F_Eq + fcopAC(iR,iZ)*dlog(fcopAC(iR,iZ))*avpola(1,iR,iZ,1)*cteAc
F_Eq = F_Eq + fcopANC(iR,iZ)*dlog(fcopANC(iR,iZ))*avpola(1,iR,iZ,1)*cteAc
F_Eq = F_Eq - fcopANC(iR,iZ)*dlog(expmuHplus)*avpola(1,iR,iZ,1)*cteAc
F_eq = F_eq + fcopANC(iR,iZ)*dlog(Ka(1))*avpola(1,iR,iZ,1)*cteAc

!eq acido base de la base

F_Eq = F_Eq + fmolC(iR,iZ)*dlog(fmolC(iR,iZ))*avpolb(1,iR,iZ,2)*cteBasic
F_Eq = F_Eq + fmolNC(iR,iZ)*dlog(fmolNC(iR,iZ))*avpolb(1,iR,iZ,2)*cteBasic
F_Eq = F_Eq - fmolNC(iR,iZ)*dlog(expmuOHmin)*avpolb(1,iR,iZ,2)*cteBasic
F_Eq = F_Eq + fmolNC(iR,iZ)*dlog(Kb(1))*avpolb(1,iR,iZ,2)*cteBasic

!eq par ionico
F_parCopA = 0.0
F_parCopA = F_parCopA - fcopAion(iR,iZ)*dlog(expmupos)*avpola(1,iR,iZ,1)*cteAc
!F_parCopA = F_parCopA + fcopAion(iR,iZ)*dlog(Kcopion/vpos)*avpola(1,iR,iZ,1)*cteAc
F_parCopA = F_parCopA - fcopAion(iR,iZ)*dlog(Kcopion)*avpola(1,iR,iZ,1)*cteAc
F_parCopA = F_parCopA + fcopAion(iR,iZ)*dlog(fcopAion(iR,iZ))*avpola(1,iR,iZ,1)*cteAc

F_eq = F_eq + F_parCopA 

F_parMol = 0.0
F_parMol = F_parMol - fmolion(iR,iZ)*dlog(expmuneg)*avpolb(1,iR,iZ,2)*cteBasic
!F_parMol = F_parMol + fmolion(iR,iZ)*dlog(Kmolion/vneg)*avpolb(1,iR,iZ,2)*cteBasic
F_parMol = F_parMol - fmolion(iR,iZ)*dlog(Kmolion)*avpolb(1,iR,iZ,2)*cteBasic
F_parMol = F_parMOl + fmolion(iR,iZ)*dlog(fmolion(iR,iZ))*avpolb(1,iR,iZ,2)*cteBasic

F_eq = F_eq + F_parMol

! eq de as CopA-mol
F_Asoc = 0.0
F_Asoc = F_Asoc - dlog(Kcopmol)*fASmol(iR,iZ)*avpola(1,iR,iZ,1)*cteAc
if(fASmol(iR,iZ).gt.1.0d-10) F_Asoc = F_Asoc + fASmol(iR,iZ)*dlog(fASmol(iR,iZ))*avpola(1,iR,iZ,1)*cteAc
if(fAScopA(iR,iZ).gt.1.0d-10) F_Asoc = F_Asoc + fAScopA(iR,iZ)*dlog(fAScopA(iR,iZ))*avpolb(1,iR,iZ,2)*cteBasic
!F_eq = F_eq - dlog(Kcopmol)*(fASmol(iR,iZ)*avpola(1,iR,iZ,1)*cteAc+fAScopA(iR,iZ)*avpolb(1,iR,iZ,2)*cteBasic)

if((avpola(1,iR,iZ,1).gt.1.0d-10).and.(fASmol(iR,iZ).gt.1.0d-10)) F_Asoc = F_Asoc - & 
        &fASmol(iR,iZ)*avpola(1,iR,iZ,1)*(dlog(avpola(1,iR,iZ,1)*fASmol(iR,iZ)*vcopmol*vsol) - 1.0)*cteAc

F_eq = F_eq + F_Asoc

print*, 'FMOLION', F_parMol, 'FCOPION', F_parCopA, 'FCOPMOL', F_Asoc

enddo
enddo

Free_Energy = Free_Energy + F_Eq     

!do im = 1, N_monomer                                                                       
!do iC  = 1, ncells                                               
!if(zpol(im).ne.0) then
!F_Eq = F_Eq + fdis(im, iC)*dlog(fdis(im, iC))
!& *avpol_monom(im, iC)/vpol     
!& *(dfloat(indexa(iC,1))-0.5)                                     
!F_Eq = F_Eq + (1.0-fdis(im, iC))                                     
!& *dlog(1.0-fdis(im, iC))*avpol_monom(im, iC)/vpol                      
!& *(dfloat(indexa(iC,1))-0.5)                                     

!F_Eq = F_Eq + (1.0-fdis(im, iC))*
!& dlog(K0(im))*avpol_monom(im, iC)/vpol     
!& *(dfloat(indexa(iC,1))-0.5)                                     

!select case (zpol(im))
!case (-1) ! acid
!F_Eq = F_Eq + (1.0-fdis(im, iC))                                     
!& *(-dlog(expmuHplus))*avpol_monom(im, iC)/vpol                     
!& *(dfloat(indexa(iC,1))-0.5)
!case (1) ! base
!F_Eq = F_Eq + (1.0-fdis(im, iC))                                     
!& *(-dlog(expmuOHmin))*avpol_monom(im, iC)/vpol                     
!& *(dfloat(indexa(iC,1))-0.5)
!end select
!endif
!enddo                                                            
!enddo                                                            
!F_eq = F_eq *delta**3/vsol*2*pi                                  
!Free_Energy = Free_Energy + F_Eq                                 

! 8.vdW ! Ojo, los  son negativos => atraccion                         

F_VdW = 0.0                                                      

do MC = 1, Ncomp
do NC = 1, Ncomp
do iR = 1, dimR
do iZ = 1, dimZ
  do is = 1, Npoorsv
  do js = 1, Npoorsv
    do jR= 1, dimR
    do jZ= -Xulimit, Xulimit    
kZ = jZ + iZ
if(PBCflag.eq.1)kkZ = PBCSYMI (kZ,dimZ)
if(PBCflag.eq.2)kkZ = PBCREFI (kZ,dimZ)
F_vdW (is,js) = F_vdW(is,js) &
        - 0.5*Xu(iR,jR,jZ,is,js)*avpol(is,iR,iZ,NC)*avpol(js,jR,kkZ,MC)  &
        *st(is,js)/(vpol(is)*vpol(js)*vsol**2)*jacobian(iR)*deltaR*deltaZ
    enddo ! jZ
    enddo ! jR
  enddo ! js
  enddo ! is
enddo ! iZ          
enddo ! iR                                             
enddo ! NC
enddo ! MC


do is=1,Npoorsv
  do js=1,Npoorsv
    Free_Energy = Free_Energy + F_vdW (is,js)
  enddo
enddo                                

!! 9. Electrostati -- no charge on surfaces                            
!LOKE

F_electro = 0.0                                                  

 
do iR = 1, dimR
do iZ = 1, dimZ

  jZp= iZ+1
  jZm= iZ-1

if(PBCflag.eq.1) then
 iZp= PBCSYMI(jZp,dimZ)
 iZm= PBCSYMI(jZm,dimZ)
else if(PBCflag.eq.2) then
 iZp= PBCREFI(jZp,dimZ)
 iZm= PBCREFI(jZm,dimZ)
endif

  gradphi2=((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2 + ((phi(iR,iZp)-phi(iR,iZm))/2.0/deltaZ)**2
  F_electro = F_electro + (xcharge(iR,iZ)*phi(iR,iZ) - wperm*epsfcn(iR,iZ)/2.0*gradphi2)*jacobian(iR)*deltaR*deltaZ
enddo
enddo

Free_Energy = Free_Energy + F_electro                            

if(rank.eq.0)print*, 'Free Energy, method I : ', Free_Energy

! Method II                                                          

Free_Energy2 = 0.0                                               

sumpi = 0.0                                                   
sumrho=0.0                                                    
sumel=0.0                                                     
sumdiel=0.0                                                     
sumpol = 0.0
sumeq = 0.0

do iR= 1, dimR                                                
do iZ= 1, dimZ
  sumpi = sumpi+dlog(xsol(iR,iZ))*jacobian(iR)                      
  sumpi = sumpi-dlog(xsolbulk)*jacobian(iR)

  sumrho = sumrho + (-xsol(iR,iZ)*jacobian(iR)) ! sum over  rho_i i=+,-,si
  sumrho = sumrho - (-xsolbulk)*jacobian(iR) ! sum over  rho_i i=+,-,si
  sumrho = sumrho + (-avpos(iR,iZ)/vpos)*jacobian(iR)
  sumrho = sumrho - (-xposbulk/vpos)*jacobian(iR)
  sumrho = sumrho + (-avneg(iR,iZ)/vneg)*jacobian(iR)
  sumrho = sumrho - (-xnegbulk/vneg)*jacobian(iR)
  sumrho = sumrho + (-avHplus(iR,iZ)*jacobian(iR))
  sumrho = sumrho - (-xHplusbulk*jacobian(iR))
  sumrho = sumrho + (-avOHmin(iR,iZ)*jacobian(iR))
  sumrho = sumrho - (-xOHminbulk*jacobian(iR))

enddo
enddo

do NC = 1, NComp
do iR=1,maxntotcounterR
do iZ=1,maxntotcounterZ
  sumrho = sumrho + (-xpol(iR,iZ,NC)*vsol*jacobian(iR)) ! sum over  rho_i i=+,-,si
enddo
enddo
enddo

! Chemical equilibrium LEO
do iR=1,dimR
do iZ=1,dimZ
 !test = test + (avpola(is,iR,iZ,Nc)*fboundA(is,ir,iz))-(avpolb(is,iR,iZ,Nc)*fboundB(is,ir,iz))
   sumeq = sumeq + avpola(1,iR,iZ,1)*fASmol(iR,iZ)/(vpol_a(1))*jacobian(iR) ! sum over  rho_i i=+,-,si
enddo
enddo

! sumel

do iR = 1, dimR
do iZ = 1, dimZ

  jZp= iZ+1
  jZm= iZ-1

if(PBCflag.eq.1) then
 iZp= PBCSYMI(jZp,dimZ)
 iZm= PBCSYMI(jZm,dimZ)
else if(PBCflag.eq.2) then
 iZp= PBCREFI(jZp,dimZ)
 iZm= PBCREFI(jZm,dimZ)
endif

  gradphi2=((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2 + ((phi(iR,iZp)-phi(iR,iZm))/2.0/deltaZ)**2
  sumel = sumel - wperm*epsfcn(iR,iZ)/2.0*gradphi2*jacobian(iR)*deltaR*deltaZ

enddo ! iZ
enddo ! R


! contribution from dielecrtric

do iR = 1, dimR
do iZ = 1, dimZ

  jZp = iZ + 1
  jZm = iZ - 1

if(PBCflag.eq.1) then
 iZp= PBCSYMI(jZp,dimZ)
 iZm= PBCSYMI(jZm,dimZ)
else if(PBCflag.eq.2) then
 iZp= PBCREFI(jZp,dimZ)
 iZm= PBCREFI(jZm,dimZ)
endif

  gradphi2 = ((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2 + ((phi(iR,iZp)-phi(iR,iZm))/2.0/deltaZ)**2
  sumdiel = sumdiel + 0.5*wperm*dielpol(iR,iZ)*Depsfcn(iR,iZ)*gradphi2*jacobian(iR)*vsol

enddo

enddo

Free_Energy2 = (sumpi + sumrho + sumdiel + sumeq)/vsol*deltaR*deltaZ + sumel                      

do is=1,Npoorsv
  do js=1,Npoorsv
    Free_Energy2 = Free_Energy2 - F_vdW(is,js)
  enddo
enddo

do NC = 1,Ncomp
mupol = dlog(xpol(1,1,NC))-dlog(q(1,1,NC))
do iR = 1, maxntotcounterR
do iZ = 1, maxntotcounterZ
  sumpol = sumpol + xpol(iR,iZ,NC)*mupol*jacobian(iR)*deltaR*deltaZ
enddo
enddo
enddo

Free_Energy2 = Free_Energy2 + sumpol 

if(rank.eq.0)print*, 'Free Energy, method II: ', Free_Energy2

if(rank.eq.0) then                                                                 
  write(301,*) npol, Free_energy/npol                       
  write(302,*) npol, F_Mix_s/npol                           
  write(303,*) npol, F_Mix_pos/npol                       
  write(304,*) npol, F_Mix_neg/npol                         
  write(305,*) npol, F_Mix_Hplus/npol                       
  write(306,*) npol, F_Mix_OHmin/npol                       
  write(307,*) npol, F_Conf/npol                            
  write(308,*) npol, F_Eq/npol                            
! write(313,*)counter, counter2, F_Eq_P                              
  do is=1,Npoorsv
    do js=1,Npoorsv
       write(10000*is+js,*) npol, F_vdW(is,js)/npol                             
    enddo
  enddo
  write(3110,*) npol, F_electro/npol                         
  write(312,*) npol, Free_energy2/npol                      
  write(313,*) npol, F_Mix_p/npol                          
  write(314,*) npol, F_Uchain/npol
endif                           

! Save end-to-end distances         

!ii = rank+1
!endtoend_av_tosend = 0.0
!endtoend_av_tosend(ii) = endtoend_av

!call MPI_REDUCE(endtoend_av_tosend, endtoend_av_all
!& , N_chains,
!&   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

!if(rank.eq.0) then
!open(unit=502, file='endtoend.dat')                          
!do ii = 1, N_chains
!write(502,*)ii, endtoend_av_all(ii)
!enddo                      
!endif

continue

end                                                              




