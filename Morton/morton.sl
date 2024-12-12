#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:01:59
#SBATCH --output=morton.out
#SBATCH -Aanakano_429

module load gcc
module load openmpi

mpirun -bind-to none -n 2 ./morton_hmd