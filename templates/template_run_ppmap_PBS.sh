#!/bin/bash
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l place=scatter:exc1
#PBS -o run_ppmap.log
#PBS -e run_ppmap.err
#PBS -N run_ppmap
#PBS -l walltime=72:00:00
#PBS -q workq
#PBS -P [JOBID]


# latest intel compilers, mkl and intel-mpi

module purge
module load 

ulimit -s unlimited
ulimit -c 0

code=
MYPATH=

NCPUS="40"

cd ${MYPATH}
start="$(date +%s)"
  for n in $NCPUS; do
      echo Running PPMAP with OMP_NUM_THREADS=$n 
      export OMP_NUM_THREADS=$n
      ${code}   <field name>   <first field #>   <last field #> 
      echo PPMAP finished
  done
#
