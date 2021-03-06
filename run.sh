#!/bin/bash

export LD_LIBRARY_PATH="$(nc-config --prefix)/lib:$(nc-config --prefix)/lib64:$(nf-config --prefix)/lib:$(nf-config --prefix)/lib64:$LD_LIBRARY_PATH"

NUM_THREADS=1
NUM_PROCS=8

mpirun -np $NUM_PROCS -x OMP_NUM_THREADS=$NUM_THREADS ./gp