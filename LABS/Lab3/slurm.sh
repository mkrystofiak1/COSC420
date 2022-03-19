#!/bin/bash

# File: myproject.sh

#SBATCH --job-name=GaussJordan
#SBATCH --ntasks=12
#SBATCH --nodelist=hpcl5-1,hpcl6-[1-3]
##SBATCH --ntasks-per-node=1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --output=out/%j.log

module load mpi/mpich-3.2-x86_64

valgrind mpiexec ~/COSC420/Lab3/l3
