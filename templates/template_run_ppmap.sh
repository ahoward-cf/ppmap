#!/bin/bash --login
#SBATCH --job-name=run_ppmap
#SBATCH -o run_ppmap.log
#SBATCH -e run_ppmap.err
#SBATCH -t 3-00:00
#SBATCH -p compute
#SBATCH -n 40
#SBATCH --exclusive

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
