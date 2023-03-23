module transgauche
real*8, allocatable :: trans(:,:)
integer*1, allocatable :: Ntrans(:,:,:)
real*8, allocatable :: Uchain(:,:)
endmodule


module mkinsol
double precision, allocatable :: pp(:)
endmodule


module mkai
integer, allocatable :: segpoorsv(:,:)
integer Npoorsv ! number of different types of poor solvent
REAL*8, allocatable ::  xtotal(:,:,:)
real*8, allocatable :: st(:,:)
real*8, allocatable :: dimf(:,:), dimfkais(:,:)
real*8, allocatable :: Xu(:,:,:,:,:)
integer Xulimit, Xulimitkais
integer MCfactor, MCfactorkais
integer flagkai, curvkais
real*8 dimRkais
endmodule

module globals
integer vtkflag, maxT
integer PBCflag
real*8, allocatable :: xflag(:)
real*8 Na
parameter (Na=6.02d23)
integer Ncomp
integer, parameter :: base = 80
real*8 lseg
real*8 lsegkai
INTEGER actionflag
real*8 npolini, npolfirst, npollast, npolstep
real*8, allocatable :: npolratio(:) 
real*8 npol
real*8 error              ! error imposed accuaracy
real*8 infile             ! inputfile control variable for reading input files  value 0,1
CHARACTER nada
real*8 norma
INTEGER adsmax
integer ntot, dimRini, dimR, dimZ, maxntotR, maxntotZ, maxntotcounterR, maxntotcounterZ, minntotR, minntotZ !lattice sites
integer Rini_kais, Rfin_kais, minntotRkais, maxntotRkais ! kai limits
real*8, allocatable :: avpol(:,:,:,:) ! volume fraction of chains 
real*8, allocatable :: avpola(:,:,:,:), avpolb(:,:,:,:) ! volume fraction of acid and basic segments 
real*8, allocatable :: xpol(:,:,:) ! volume fraction polymers already adsorbed
real*8, allocatable :: xsol(:,:)
INTEGER totalcuantas, cuantas
integer curvature
real*8, allocatable :: Ug(:), Ut(:)

real*8, allocatable :: epsfcn(:,:)
real*8, allocatable :: dielpol(:,:)
real*8, allocatable :: Depsfcn(:,:)

integer first, last

integer*2, allocatable :: innR(:,:,:,:), innZ(:,:,:)

integer iter              ! counts number of iterations

! integer, parameter :: ncha_max = 700
integer ncha_max

endmodule

module partfunc
real*8, allocatable :: q(:,:,:)
real*8, allocatable ::  sumprolnpro(:,:,:), sumprouchain(:,:,:)
endmodule

module layer
real*8 deltaR 
real*8 deltaZ 
endmodule

module volume
real*8, allocatable :: vchain(:)
real*8 vsol, vpos, vneg, r_neg, r_pos
real*8, allocatable :: vpol(:), vpol_a(:), vpol_b(:)
endmodule

module mcharge
real*8 dielP
integer electroflag
integer Nacids, Nbasics
integer, allocatable :: acidtype(:,:), basictype(:,:)
real*8, allocatable :: pKa(:), pKb(:), Ka(:), Kb(:)
real*8, allocatable :: phi(:,:), avpos(:,:), avneg(:,:), avHplus(:,:), avOHmin(:,:), xcharge(:,:), fAmin(:,:,:), fBHplus(:,:,:)
real*8 Csalt, wperm, rhosalt, expmupos, expmuneg, xposbulk, xnegbulk
real*8 cHplus, cOHmin, pHbulk, pOHbulk, pKw, xHplusbulk, xOHminbulk, expmuHplus, expmuOHmin
endmodule


module bulk
REAL*8 expmupol
real*8 xsolbulk, phibulkpol           ! volume fraction of solvent in bulk
endmodule

module seed1
integer seed, last_seed              ! seed for random number generator
endmodule


module longs
integer entflag
integer, allocatable :: long(:)            ! length of polymer
integer maxlong
integer, allocatable :: nbranches(:)
integer maxnbranches
integer, allocatable :: long_branches(:)
integer, allocatable :: branch_pos(:,:), branch_long(:,:) ! position and lenght of branches read from input
integer, allocatable :: torsionstate(:,:)
endmodule

module pis
real*8 pi
endmodule

module matrices
real*8 tt(3,3),tp(3,3),tm(3,3)
endmodule

module senos
real*8 ta,sitheta,cotheta,siphip,cophip
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
integer status(MPI_STATUS_SIZE)
endmodule

