#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:01:59
#SBATCH --output=test_pmd.out
#SBATCH -Aanakano_429

module load openmpi
module load gcc

mpirun -bind-to none -n 2 ./test_pmd

