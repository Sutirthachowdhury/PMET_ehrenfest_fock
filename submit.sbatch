#!/bin/bash
#SBATCH -p action -A action
#SBATCH -o rpmd.out
#SBATCH --mem-per-cpu=1GB
#SBATCH -t 120:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#module load lapack/3.7.0/b1

time ./rpmd.exe
