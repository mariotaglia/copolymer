
module mkinsol
double precision, allocatable :: pp(:)
endmodule

module globals

integer, parameter :: base = 80

real*8 lseg
real*8, allocatable :: Xu(:,:)
integer Xulimit
REAL*8 sts(1000)
INTEGER nst, actionflag

real*8 npolini, npolfirst, npollast, npolstep 

real*8 npol
real*8 error              ! error imposed accuaracy
real*8 infile             ! inputfile control variable for reading input files  value 0,1
CHARACTER nada

real*8 norma
INTEGER adsmax

integer ntot, maxntot, maxntotcounter_ini, maxntotcounter ! lattice sites
real*8, allocatable :: avpol(:,:) ! volume fraction polymers already adsorbed
real*8, allocatable :: xpol(:) ! volume fraction polymers already adsorbed
real*8, allocatable :: xsol(:)
INTEGER cuantas
integer curvature

integer*2, allocatable :: in1n(:,:,:)
integer*2, allocatable :: in2n(:,:,:)
integer, allocatable ::  maxpos(:,:)
integer, allocatable ::  minpos(:,:)


real*8, allocatable :: eps(:)
real*8 eps1

integer iter              ! counts number of iterations

REAL*8, allocatable ::  xtotal(:)
real*8 st
integer, parameter :: ncha_max = 700
endmodule

module partfunc
real*8, allocatable :: q(:)
real*8 qall
real*8 sumprolnproall
endmodule

module layer
real*8, parameter :: delta = 0.2
endmodule

module volume
real*8 vpol, vsol
endmodule


module bulk
REAL*8 expmupol
real*8 xsolbulk, phibulkpol           ! volume fraction of solvent in bulk
endmodule

module seed1
integer seed              ! seed for random number generator
endmodule


module longs
integer long, long1            ! length of polymer
endmodule

module pis
real*8 pi
endmodule

module matrices
real*8 tt(3,3),tp(3,3),tm(3,3)
endmodule

module senos
real*8 sitheta,cotheta,siphip,cophip
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
endmodule

