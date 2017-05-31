subroutine calc_free_energy(counter, counter2)

use globals
use partfunc
use pis
use MPI
use bulk
use layer
use volume
implicit none                                 

double precision Factorcurv                   

real*8 Free_energy, F_Mix_s, F_Mix_pos, F_mix_p        
real*8 F_Mix_neg, F_Mix_Hplus                 
real*8 Free_energy2, sumpi, sumrho, sumel, sum, sumpol, mupol
real*8 F_Mix_OHmin, F_Conf       
real*8  F_Conf2, F_Conf_temp2, F_Eq, F_Eq_P, F_vdW, F_eps, F_electro                            

integer counter, counter2                                    

integer iC, jC                                  

real*8, external :: jacobian

if(rank.eq.0) then 
print*, 'Starting free energy calculation...'
! open files
open(unit=301, file='F_tot.dat')                           
open(unit=302, file='F_mixs.dat')                          
!open(unit=303, file='F_mixpos.dat')                        
!open(unit=304, file='F_mixneg.dat')                        
!open(unit=305, file='F_mixH.dat')                          
!open(unit=306, file='F_mixOH.dat')                         
open(unit=307, file='F_conf.dat')                          
!open(unit=308, file='F_eq.dat')                            
!open(unit=313 ,file='F_eq_P.dat')                            
open(unit=309, file='F_vdW.dat')                           
!open(unit=310, file='F_eps.dat')                           
!open(unit=311, file='F_electro.dat')                       
open(unit=312, file='F_tot2.dat')                          
open(unit=313, file='F_mixp.dat')                          
endif


!---------------------------------------------------------------------
!        Calculate Free Energy                                        
!-------------------------------------------------------------------- 

Free_Energy = 0.0                                                
Free_Energy2 = 0.0                                               

! 1. Solvent translational entropy
F_Mix_s = 0.0                                                    

do iC = 1, ntot                                                
F_Mix_s = F_Mix_s + xsol(iC)*(dlog(xsol(iC))-1.0)*jacobian(iC)*delta/vsol                                      
F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)*jacobian(iC)*delta/vsol                              
enddo                                                            
Free_Energy = Free_Energy + F_Mix_s                              


! 1. Solvent translational entropy
F_Mix_p = 0.0                                                    

do iC = 1, maxntot                                                
F_Mix_p = F_Mix_p + xpol(iC)*(dlog(xpol(iC))-1.0)*jacobian(iC)*delta                                      
enddo                                                            
Free_Energy = Free_Energy + F_Mix_p            

! 2. Pos ion translational entropy

!F_Mix_pos = 0.0                                                  

!do iC = 1, ntot                                               
!F_Mix_pos = F_Mix_pos + xpos(iC)*(dlog(xpos(iC)/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*jacobian(iC)*delta/vsol
!F_Mix_pos = F_Mix_pos - xposbulk*(dlog(xposbulk/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))*jacobian(iC)*delta/vsol                               
!enddo                                                            
!Free_Energy = Free_Energy + F_Mix_pos                            

!3. Neg ion translational entropy

!F_Mix_neg = 0.0                                                  

!do iC = 1, ntot                                                
!F_Mix_neg = F_Mix_neg + xneg(iC)*(dlog(xneg(iC)/vsalt)-1.0-dlog(expmuneg) + dlog(vsalt))*jacobian(iC)*delta/vsol       
!F_Mix_neg = F_Mix_neg - xnegbulk*(dlog(xnegbulk/vsalt)-1.0-dlog(expmuneg) + dlog(vsalt))*jacobian(iC)*delta/vsol       
!enddo                                                            
!Free_Energy = Free_Energy + F_Mix_neg                            


! 4. H+ translational entropy

!F_Mix_Hplus = 0.0                                                

!do iC = 1, ntot                                      
!F_Mix_Hplus = F_Mix_Hplus + xHplus(iC)*(dlog(xHplus(iC))-1.0-dlog(expmuHplus))*jacobian(iC)*delta/vsol                                             
!F_Mix_Hplus = F_Mix_Hplus - xHplusbulk*(dlog(xHplusbulk)-1.0-dlog(expmuHplus))*jacobian(iC)*delta/vsol                                             
!enddo                                                            
!Free_Energy = Free_Energy + F_Mix_Hplus                          

! 5. OH- translational entropy 

!F_Mix_OHmin = 0.0                                                

!do iC = 1, ntot                                       
!F_Mix_OHmin = F_Mix_OHmin + xOHmin(iC)*(dlog(xOHmin(iC))-1.0-dlog(expmuOHmin))*jacobian(iC)*delta/vsol           
!F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))*jacobian(iC)*delta/vsol                                             
!enddo                                                            
!Free_Energy = Free_Energy + F_Mix_OHmin                          

! 6. Polymer conformational entropy                                         


F_conf = 0 
do iC = 1, maxntot 
F_Conf = F_conf + (sumprolnpro(iC)/q(iC)-dlog(q(iC)))*jacobian(iC)*delta*xpol(iC)
enddo 
Free_Energy = Free_Energy + F_Conf

! 7. Chemical Equilibrium                                              

!F_Eq = 0.0                                                       
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

do iC = 1, ntot
do jC = 1, ntot                                         

F_vdW = F_vdW - 0.5*Xu(iC,jC)*avpol(iC,2)*avpol(jC,2)*st/((vpol*vsol)**2)*jacobian(iC)*delta

!F_vdW = F_vdW - 0.5000*delta**3*xtotal2(ii,iC)*      
!&       xtotal2(iii,Xulist_cell(iC, iiC))*                    
!&       Xulist_value(iC,iiC)*st_matrix(ii, iii)*st
!&       /((vpol*vsol)**2)                   
!&       *(dfloat(indexa(iC,1))-0.5)*2*pi               

enddo ! iC          
enddo ! jC                                             

Free_Energy = Free_Energy + F_vdW                                

!! 9. Electrostati -- no charge on surfaces                            
!F_electro = 0.0                                                  
!do iC  = 1, ncells                                               
!F_electro = F_electro + delta**3*psi2(iC)*qtot(iC)/2.0/vsol      
!&               *(dfloat(indexa(iC,1))-0.5)*2*pi                  
!enddo                                                            
!Free_Energy = Free_Energy + F_electro                            

if(rank.eq.0)print*, 'Free Energy, method I: ', Free_Energy

! Method II                                                          

Free_Energy2 = 0.0                                               

sumpi = 0.0                                                   
sumrho=0.0                                                    
sumel=0.0                                                     
sumpol = 0.0

do iC=1,ntot                                                
sumpi = sumpi+dlog(xsol(iC))*jacobian(iC)                      
sumpi = sumpi-dlog(xsolbulk)*jacobian(iC)

sumrho = sumrho + (-xsol(iC)*jacobian(iC)) ! sum over  rho_i i=+,-,si
sumrho = sumrho - (-xsolbulk)*jacobian(iC) ! sum over  rho_i i=+,-,si
enddo


do iC=1,maxntot                                                
sumrho = sumrho + (-xpol(iC)*vsol*jacobian(iC)) ! sum over  rho_i i=+,-,si
enddo

!do iC=1,ntot                                                
!sumel = sumel - qtot(iC)*psi2(iC)/2.0
!&               *(dfloat(indexa(iC,1))-0.5)*2*pi
!sumel = sumel + proteinqC(iC)*psi2(iC)*vsol                   
!&               *(dfloat(indexa(iC,1))-0.5)*2*pi   
!enddo                                                            

Free_Energy2 = (sumpi + sumrho + sumel)/vsol*delta                               
Free_Energy2 = Free_Energy2 - F_vdW

mupol = dlog(xpol(1))-dlog(q(1))
do iC = 1, maxntot
sumpol = sumpol + xpol(iC)*mupol*jacobian(iC)*delta
enddo
Free_Energy2 = Free_Energy2 + sumpol

if(rank.eq.0)print*, 'Free Energy, method II: ', Free_Energy2

if(rank.eq.0) then                                                                 
write(301,*)counter, counter2, Free_energy/npol                       
write(302,*)counter, counter2, F_Mix_s/npol                           
!write(303,*)counter, counter2, F_Mix_pos                         
!write(304,*)counter, counter2, F_Mix_neg                         
!write(305,*)counter, counter2, F_Mix_Hplus                       
!write(306,*)counter, counter2, F_Mix_OHmin                       
write(307,*)counter, counter2, F_Conf/npol                            
!write(308,*)counter, counter2, F_Eq                              
!write(313,*)counter, counter2, F_Eq_P                              
write(309,*)counter, counter2, F_vdW/npol                             
!write(310,*)counter, counter2, F_eps                          
!write(311,*)counter, counter2, F_electro                         
write(312,*)counter, counter2, Free_energy2/npol                      
write(313,*)counter, counter2, F_Mix_p/npol                          
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




