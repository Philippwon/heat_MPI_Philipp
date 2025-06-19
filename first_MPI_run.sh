#!/usr/bin/env bash

# Job Name and Files (also --job-name)
#SBATCH -J heat
#SBATCH --mail-user=ge57jij@mytum.de
#Output and error (also --output, --error):
#SBATCH -o job2.out
#SBATCH -e job2.out
# Wall clock limit:
#SBATCH --time=00:03:00
#SBATCH --account=h039v
#SBATCH --partition=test
#SBATCH --nodes=15
#--------------------------------------
module load slurm_setup



# mpirun -np 21 ./heat jacobi.dat heat_neu.ppm 3 7
mpiexec -n 15 ./heat_MPI jacobi.dat heat_neu.ppm 5 3
