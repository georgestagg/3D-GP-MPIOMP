#!/bin/bash
#SBATCH --ntasks=128
#SBATCH -c 1

module load netCDF-Fortran/4.4.4-goolf-1.5.14
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="$(nc-config --prefix)/lib:$(nc-config --prefix)/lib64:$(nf-config --prefix)/lib:$(nf-config --prefix)/lib64:$LD_LIBRARY_PATH"

mpirun ./gp
