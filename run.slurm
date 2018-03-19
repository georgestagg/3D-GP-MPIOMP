#!/bin/bash
#SBATCH --ntasks=128
#SBATCH -c 1

module use /nobackup/proj/mbsmsp/modules
module load netcdf-fortran/4.4.4
module load FFTW/3.3.6-gompi-2017b
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="$(nc-config --prefix)/lib:$(nc-config --prefix)/lib64:$(nf-config --prefix)/lib:$(nf-config --prefix)/lib64:$LD_LIBRARY_PATH"

mpirun -x OMP_NUM_THREADS=1 ./gp
