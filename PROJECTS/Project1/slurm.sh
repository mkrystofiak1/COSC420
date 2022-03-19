#!/bin/bash

# File: myproject.sh

#SBATCH --job-name=democluster
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=12
##SBATCH --ntasks=50
#SBATCH --mem=4gb
#SBATCH --time=24:00:00
#SBATCH --output=out/cracker.log

module load mpi/mpich-3.2-x86_64

mpiexec ~/COSC420/Projects/Project1/cracker
