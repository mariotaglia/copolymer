
module mkinsol
double precision, allocatable :: pp(:)
endmodule

module global
real*8 lseg
real*8, allocatable :: Xu(:,:)
real*8 sumXu11, sumXu12
integer Xulimit
REAL*8 sts(100), npols(100)
real*8 npol
INTEGER nst, nnpol
real*8 error              ! error imposed accuaracy
real*8 infile             ! inputfile control variable for reading input files  value 0,1
CHARACTER nada

real*8 norma
integer cuantas          ! number of polymer configuration or  bound sequences

integer ntot ! lattice sites
real*8, allocatable :: avpol(:,:) ! volume fraction polymers already adsorbed

integer*4, allocatable :: in1n(:,:,:)
integer*4, allocatable :: in2n(:,:,:)

real*8 sigma

real*8, allocatable :: eps(:)
real*8 eps1

integer iter              ! counts number of iterations

REAL*8, allocatable ::  xtotal(:)
real*8 st
endmodule

module partfunc
real*8 q
endmodule

module layer
real*8, parameter :: delta = 0.5
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
integer long(2)            ! length of polymer

INTEGER long1,long2,maxlong
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

module posmk
real*8, allocatable :: current(:,:)
integer*2, allocatable :: nextbead(:)
integer*2, allocatable :: firstcell(:,:,:)
integer, parameter :: mcube = 100
integer, parameter :: calq = 0
real*8, parameter :: qprob0 = 0.6933
integer, parameter :: nearbonds = 5
endmodule posmk

