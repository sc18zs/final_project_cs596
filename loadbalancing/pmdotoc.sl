#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1  # Run 4 tasks per node (using all CPUs)
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00  # Set a 10-minute time limit to prevent it from running too long
#SBATCH --output=newout.out
#SBATCH --error=pmdotoc.err  # Separate error log file

mpirun -bind-to none -n 2 ./pmdotoc 