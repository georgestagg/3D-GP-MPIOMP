#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/.fs/netcdf/lib:/data/.fs/netcdf-fortran/lib

export OMP_NUM_THREADS=2
export MPI_NUM_NODES=4

mpirun -np $MPI_NUM_NODES ./gp
