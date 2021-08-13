#!/bin/bash

#SBATCH --mail-user=abento@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=0-15:0:00
#SBATCH --mem=58gb
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=Omega_RA_0.1

module load r

for RA in 0.1
do
    for gridsize in 21
    do
      label="Omega_RA_${RA}"
      nohup Rscript epicon_script.R 96 $gridsize 1 "nosocial" "decentralized" 5.166426 7.513051 3.548544 2.6 5.1 58000 $RA 1 "eqm_$label" "R0" 1 0.8554632 1 1 "-672827.916835845" "Omega" > logs/eqm_$label.log 2>&1

      nohup Rscript epicon_script.R 96 $gridsize 1 "nosocial" "planner" 5.166426 7.513051 3.548544 2.6 5.1 58000 $RA 1 "planner_$label" "R0" 2 0.8554632 1 1 "-672827.916835845" "Omega" > logs/planner_$label.log 2>&1
    done
done 
