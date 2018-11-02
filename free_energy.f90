subroutine calc_free_energy(counter, counter2)

use globals
use mcharge
use partfunc
use pis
use MPI
use bulk
use layer
use volume
use mkai
use longs
implicit none                                 

double precision Factorcurv                   

real*8 gradphi2
real*8 Free_energy, F_Mix_s, F_Mix_pos, F_mix_p        
real*8 F_Mix_neg, F_Mix_Hplus                 
real*8 Free_energy2, sumpi, sumrho, sumel, sum, sumpol, mupol, sumdiel
real*8 F_Mix_OHmin, F_Conf, F_Uchain     
real*8  F_Conf2, F_Conf_temp2, F_Eq, F_Eq_P, F_vdW(Npoorsv,Npoorsv), F_electro                            
!real*8 F_mup
integer counter, counter2                                    
real*8 gradphi2
integer i, iC, jC, is, js, iR, iZ, jR, jZ, kZ, kkZ, jZp, jZm, iZp, iZm
integer, external :: PBCSYMI
real*8, external :: jacobian

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
!open(unit=313 ,file='F_eq_P.dat')                            
do is=1,Npoorsv
do js=1,Npoorsv
write(F_vdWfilename(is,js),'(A6,BZ,I3.3,A1,I3.3,A4)')'F_vdW.',is,'.',js,'.dat'
open(unit=10000*is+js, file=F_vdWfilename(is,js) )                           
enddo
enddo
open(unit=311, file='F_electro.dat')                       
open(unit=312, file='F_tot2.dat')                          
open(unit=313, file='F_mixp.dat')                          
open(unit=314, file='F_Uchain.dat')
endif


!---------------------------------------------------------------------
!        Calculate Free Energy                                        
!-------------------------------------------------------------------- 

Free_Energy = 0.0                                                
Free_Energy2 = 0.0                                               

mupol = dlog(xpol(1))-dlog(q(1))

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

do iR = 1, maxntotcounterR                                                
do iZ = 1, maxntotcounterZ
   F_Mix_p = F_Mix_p + xpol(iR,iZ)*(dlog(xpol(iR,iZ))-1.0)*jacobian(iR)*deltaR*deltaZ ! mix entropy of chains with respect to bulk (xpolbulk=0)                                     
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
  F_Mix_OHmin = F_Mix_OHmin + avOHmin(iC)/vsol*(dlog(avOHmin(iR,iZ))-1.0-dlog(expmuOHmin))*jacobian(iR)*deltaR*deltaZ ! mix entropy and N*Mu term of OHmin          
  F_Mix_OHmin = F_Mix_OHmin - xOHminbulk/vsol*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*jacobian(iR)*deltaR*deltaZ ! with respect to bulk       
enddo                                                            
enddo

Free_Energy = Free_Energy + F_Mix_OHmin                          

! 6. Polymer conformational entropy                                         


F_conf = 0 

do iR = 1, maxntotcounterR
do iZ = 1, maxntotcounterZ
  F_Conf = F_conf + (sumprolnpro(iR,iZ)/q(iR,iZ)-dlog(q(iR,iZ)))*jacobian(iR)*deltaR*deltaZ*xpol(iR,iZ)
enddo 
enddo

Free_Energy = Free_Energy + F_Conf

F_Uchain = 0.0

do iR=1, maxntotcounterR
do iZ=1, maxntotcounterZ
  F_Uchain = F_Uchain + deltaR*deltaZ*xpol(iR,iZ)*jacobian(iR)*(sumprouchain(iR,iZ)/q(iR,iZ))
enddo
enddo

Free_Energy = Free_Energy + F_Uchain

! 7. Chemical Equilibrium                                              

F_Eq = 0.0                                                       

do iR=1,dimR
do iZ=1,dimZ

  do is=1, Nacids
    F_Eq = F_Eq + fAmin(is,iR,iZ)*dlog(fAmin(is,iR,iZ))*avpola(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_a(is)*vsol)
    F_Eq = F_Eq + (1.0-fAmin(is,iR,iZ))*dlog(1.0-fAmin(is,iR,iZ))*avpola(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_a(is)*vsol)                                     
    F_Eq = F_Eq - (1.0-fAmin(is,iR,iZ))*dlog(expmuHplus)*avpola(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_a(is)*vsol)
    F_Eq = F_Eq + (1.0-fAmin(is,iR,iZ))*dlog(Ka(is))*avpola(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_a(is)*vsol)
  enddo
  
  do is=1, Nbasics
    F_Eq = F_Eq + fBHplus(is,iR,iZ)*dlog(fBHplus(is,iR,iZ))*avpolb(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_b(is)*vsol)
    F_Eq = F_Eq + (1.0-fBHplus(is,iR,iZ))*dlog(1.0-fBHplus(is,iR,iZ))*avpolb(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_b(is)*vsol)
    F_Eq = F_Eq - (1.0-fBHplus(is,iR,iZ))*dlog(expmuOHmin)*avpolb(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_b(is)*vsol)
    F_Eq = F_Eq + (1.0-fBHplus(is,iR,iZ))*dlog(Kb(is))*avpolb(is,iR,iZ)*jacobian(iR)*deltaR*deltaZ/(vpol_b(is)*vsol)
  enddo    

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


do iR = 1, dimR
do iZ = 1, dimZ
  do is = 1, Npoorsv
  do js = 1, Npoorsv
    do jR= 1, dimR
    do jZ= -Xulimit, Xulimit    
      kZ = jZ + iZ
      kkZ = PBCSYMI (kZ,dimZ)
      F_dW (is,js) = F_vdW(is,js) - 0.5*Xu(iR,jR,jZ,is,js) &
                                      *avpol(is,iR,iZ)*avpol(js,jR,kkZ) &
                                      *st(is,js)/(vpol(is)*vpol(js)*vsol**2) &
                                      *jacobian(iR)*delta 
  
!       F_vdW = F_vdW - 0.5000*delta**3*xtotal2(ii,iC)*      
!       &       xtotal2(iii,Xulist_cell(iC, iiC))*                    
!       &       Xulist_value(iC,iiC)*st_matrix(ii, iii)*st
!       &       /((vpol*vsol)**2)                   
!       &       *(dfloat(indexa(iC,1))-0.5)*2*pi               

    enddo ! jZ
    enddo ! jR
  enddo ! js
  enddo ! is
enddo ! iZ          
enddo ! iR                                             

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
  iZp= PBCSYMI(jZp,dimZ)
  gradphi2=((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2 + ((phi(iR,iZp)-phi(iR,iZ))/deltaZ)**2
  F_electro = F_electro + (xcharge(iR,iZ)*phi(iR,iZ) - wperm*epsfcn(iR,iZ)/2.0*gradphi2*jacobian(iR)*deltaR*deltaZ
enddo
enddo

Free_Energy = Free_Energy + F_electro                            

if(rank.eq.0)print*, 'Free Energy, method I: ', Free_Energy

! Method II                                                          

Free_Energy2 = 0.0                                               

sumpi = 0.0                                                   
sumrho=0.0                                                    
sumel=0.0                                                     
sumdiel=0.0                                                     
sumpol = 0.0

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

do iR=1,maxntotcounterR
do iZ=1,maxntotcounterZ
  sumrho = sumrho + (-xpol(iR,iZ)*vsol*jacobian(iR)) ! sum over  rho_i i=+,-,si
enddo
enddo

! sumel

do iR = 1, dimR
do iZ = 1, dimZ
  jZp = iZ + 1
  iZp = PBCSYMI(jZp,dimZ)
  gradphi2=((phi(iR+1,iZ)-phi(iR,iZ))/deltaR)**2 + ((phi(iR,iZp)-phi(iR,iZ))/deltaZ)**2
  sumel = sumel - wperm*epsfcn(iR,iZ)/2.0*gradphi2*jacobian(iR)*deltaR*deltaZ
enddo
enddo

! contribution from dielecrtric

do iR = 1, dimR
do iZ = 1, dimZ
  jZp = iZ + 1
  iZp = PBCSYMI(jZp,dimZ)
  gradphi2 = ((phi(iR+1,iZ)-phi(iR,iZ))/deltaR + (phi(iR,iZp)-phi(iR,iZ)/deltaZ))**2
  sumdiel = sumdiel + 0.5*wperm*dielpol(iR,iZ)*Depsfcn(iR,iZ)*gradphi2*jacobian(iR)*vsol
enddo
enddo

Free_Energy2 = (sumpi + sumrho + sumdiel)/vsol*delta + sumel                      

do is=1,Npoorsv
  do js=1,Npoorsv
    Free_Energy2 = Free_Energy2 - F_vdW(is,js)
  enddo
enddo

do iR = 1, maxntotcounterR
do iZ = 1, maxntotcounterZ
  sumpol = sumpol + xpol(iR,iZ)*mupol*jacobian(iR)*deltaR*deltaZ
enddo
enddo

Free_Energy2 = Free_Energy2 + sumpol 
!Free_Energy2 = Free_Energy2 + F_mup

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
  write(311,*) npol, F_electro/npol                         
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

1515    continue

end                                                              




