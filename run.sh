#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/.fs/netcdf/lib:/data/.fs/netcdf-fortran/lib

export OMP_NUM_THREADS=1
export MPI_NUM_NODES=6

mpirun -np $MPI_NUM_NODES ./gp
