#!/bin/bash

# File: myproject.sh

#SBATCH --job-name=democluster
##SBATCH --ntasks=30
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH --output=out/%n.log

module load mpi/mpich-3.2-x86_64

valgrind mpiexec ~/COSC420/Lab2/main
