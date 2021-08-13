#!/bin/bash

#SBATCH --mail-user=abento@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=0-8:0:00
#SBATCH --mem=58gb
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=rhoo_8_8

module load r

for rho_o in 8
do
  for phi in 0.8
  do
    for gridsize in 10
    do
      label="rhoo_${rho_o}_phi_${phi}"
      # nohup Rscript epicon_script.R 96 $gridsize 0 "nosocial" "decentralized" 5.166426 7.513051 $rho_o 2.6 5.1 58000 0.1 1 "ngm_eqm_$label" "R0" 1 $phi > logs/eqm_$label.log 2>&1

      nohup Rscript epicon_script.R 96 $gridsize 1 "nosocial" "planner" 5.166426 7.513051 $rho_o 2.6 5.1 58000 0.1 1 "ngm_planner_$label" "R0" 2 $phi > logs/planner_$label.log 2>&1
    done
  done
done

######  Job commands go below this line #####
