#!/bin/bash
#SBATCH --ntasks=64
#SBATCH -c 1

module load netcdf-fortran
module load fftw
export LD_LIBRARY_PATH="$(nc-config --prefix)/lib:$(nc-config --prefix)/lib64:$(nf-config --prefix)/lib:$(nf-config --prefix)/lib64:$LD_LIBRARY_PATH"

srun ./gp
