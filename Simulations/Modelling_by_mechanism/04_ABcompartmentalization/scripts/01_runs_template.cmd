#!/bin/bash

#SBATCH --job-name rep_XXXreplicaXXX
#SBATCH -D XXXdirXXX
#SBATCH -n 1                   # Number of cores. For now 56 is the number max of core available
#SBATCH --mem 15Gb
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o replica_XXXreplicaXXX.out # File to which STDOUT will be written
#SBATCH -e replica_XXXreplicaXXX.out # File to which STDERR will be written 

cd XXXdirXXX

replica=XXXreplicaXXX
condition=XXXconditionXXX

python ./01_runs_${condition}_replica_${replica}.py &> replica_XXXreplicaXXX.out
rm -fvr minimization_*.XYZ
