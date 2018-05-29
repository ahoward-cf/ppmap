#!/bin/bash
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l place=scatter:excl
#PBS -o run_ppmap.log
#PBS -e run_ppmap.err
#PBS -N run_ppmap
#PBS -l walltime=3:00:00
#PBS -q workq
#PBS -P <User ID code>

# latest intel compilers, mkl and intel-mpi

module purge
module load intel

ulimit -s unlimited
ulimit -c 0

code=${HOME}/ppmap_2017-08-11/ppmap
MYPATH=/scratch/<username>

NCPUS="16"

cd ${MYPATH}
start="$(date +%s)"
  for n in $NCPUS; do
      echo Running PPMAP with OMP_NUM_THREADS=$n 
      export OMP_NUM_THREADS=$n
      ${code}   <field name>   <first field #>   <last field #> 
      echo PPMAP finished
  done
#
