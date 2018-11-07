subroutine initmpi
use MPI
use seed1
implicit none

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

if(rank.eq.0) then
  seed=435
  print*, 'I am', rank, ' and my seed is', seed
else 
  seed=0
endif

if(rank.eq.0)print*, 'Program COPOLYMER'
if(rank.eq.0)print*, 'GIT Version: ', _VERSION
end subroutine

