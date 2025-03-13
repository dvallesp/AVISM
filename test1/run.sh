#!/bin/bash


## export OMP_NUM_THREADS=16
## export OMP_STACKSIZE=2000m
## export OMP_PROC_BIND=true
## export OMP_WAIT_POLICY=active
## export KMP_BLOCKTIME=1000000
## export GOMP_CPU_AFFINITY="8-23"

export memoryuse=256000m
export vmemoryuse=256000m
export stacksize=256000m
export coredumpsize=1m

ulimit -s unlimited

./mock_voids.x ##> voids.out &

