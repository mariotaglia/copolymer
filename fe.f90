subroutine calc_free_energy(counter, counter2)

use ***

implicit none                                 
include 'mpif.h'                              
include 'MPI.h'                               

double precision Factorcurv                   

real*8 Free_energy2, sumpi, sumrho, sum, mupol                                            
real*8 counter, counter2                                    

integer iC, ii, iiC, i, jj, im, iii                                  


! vdW ! Ojo, los  son negativos => atraccion                         

F_VdW = 0.0                                                      

do i = 1, ncells                                               
do ii = 1, ncells ! loop over kai neighbors     

F_vdW = F_vdW - 0.5000*delta**3*xtotal(i)*      
&       xtotal(ii)*                    
&       Xu(i,ii)*st CHECK
&       /((vpol*vsol)**2)                   
&       *(dfloat(indexa(iC,1))-0.5)*2*pi               

enddo  ! ii                                            
enddo ! i         


c! Method II                                                          

Free_Energy2 = 0.0                                               

sumpi = 0.0                                                   
sumrho=0.0                                                    
sumel=0.0                                                     
sumfcargo = 0.0                                                                       

do iC=1,ncells                                                

sumpi =                                                    
& sumpi+dlog(xsol(iC))                         
& *(dfloat(indexa(iC,1))-0.5)*2*pi                                


sumpi = sumpi-dlog(xsolbulk)*(dfloat(indexa(iC,1))-0.5)*2*pi

sumrho = sumrho + ( - xsol(iC) - npol(iC))*(dfloat(indexa(iC,1))-0.5)*2*pi                     

sumrho = sumrho - (- xsolbulk)*(dfloat(indexa(iC,1))-0.5)*2*pi                     

sumpi = (delta**3/vsol)*sumpi                                    
sumrho = (delta**3/vsol)*sumrho                                  

Free_Energy2 = sumpi + sumrho                                      
Free_Energy2 = Free_Energy2 - F_vdW

do ii = 1, N_chains                                                
Free_Energy2 = Free_Energy2 -
& chainsperdelta(ii)*dlog(q0(ii)/shift)
enddo

if(rank.eq.0)print*, 'Free Energy, method II: ', Free_Energy2

if(rank.eq.0) then                                                                 
write(312,*)counter, counter2, Free_energy2                      
endif                           

c! Save end-to-end distances         

ii = rank+1
endtoend_av_tosend = 0.0
endtoend_av_tosend(ii) = endtoend_av

call MPI_REDUCE(endtoend_av_tosend, endtoend_av_all
& , N_chains,
&   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

if(rank.eq.0) then
open(unit=502, file='endtoend.dat')                          
do ii = 1, N_chains
write(502,*)ii, endtoend_av_all(ii)
enddo                      
endif

1515    continue

end                                                              




