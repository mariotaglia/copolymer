
module mkinsol
double precision, allocatable :: pp(:)
endmodule


module mkai
integer Npoorsv ! number of different types of poor solvent
REAL*8, allocatable ::  xtotal(:,:)
real*8, allocatable :: st(:,:)
real*8, allocatable :: Xu(:,:)
integer Xulimit
endmodule

module globals

integer, parameter :: base = 80
real*8 lseg
INTEGER actionflag
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

integer*2, allocatable :: inn(:,:,:,:)
integer, allocatable ::  maxpos(:,:)
integer, allocatable ::  minpos(:,:)


real*8, allocatable :: eps(:)
real*8 eps1

integer iter              ! counts number of iterations

integer, parameter :: ncha_max = 700

real*8, allocatable :: Ug(:), Ut(:)
real*8, allocatable :: Uchain(:)
endmodule

module partfunc
real*8, allocatable :: q(:)
real*8, allocatable ::  sumprolnpro(:)
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
integer long            ! length of polymer
integer, allocatable :: segpoorsv(:)
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

