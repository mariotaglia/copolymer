#!/bin/bash
#SBATCH --time=12:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4096      # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=_NAME
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --ntasks=1
##SBATCH --nodes=1  

module load compiler/gcc/5.4 openmpi/2.0 python/3.5
python3 ~/develop/copolymer/scripts/auto/run.py ~/develop/copolymer/assembly
