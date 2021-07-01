#!/bin/bash
#SBATCH --ntasks=64
#SBATCH -c 1

module use /nobackup/proj/mbsmsp/modules
module load netcdf-c/4.7.0 netcdf-fortran/4.4.5 fftw/3.3.8

mpirun ./gp
