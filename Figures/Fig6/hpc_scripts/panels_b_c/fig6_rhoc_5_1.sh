#!/bin/bash

#SBATCH --mail-user=abento@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=0-24:0:00
#SBATCH --mem=58gb
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=rhoc_5_1

module load r

if [[ $# -gt 0 ]]
then
  hot_start=$1
fi

if [[ -z $1 ]]
then
  hot_start=0
fi

for rho_c in 5
do
  for phi in 0.1
  do
    for gridsize in 11
    do
      label="rhoc_${rho_c}_phi_${phi}"
		nohup Rscript epicon_script_cheby.R 48 $gridsize $hot_start "nosocial" "decentralized" $rho_c 7.513051 3.548544 2.6 5.1 58000 0.1 0 "eqm_$label" "R0" 1.01 $phi 1 1 "1e07" "VSL" > logs/eqm_$label.log 2>&1

		nohup Rscript epicon_script_cheby.R 48 $gridsize $hot_start "nosocial" "planner" $rho_c 7.513051 3.548544 2.6 5.1 58000 0.1 0 "planner_$label" "R0" 1.01 $phi 1 1 "1e07" "VSL" > logs/planner_$label.log 2>&1
    done
  done
done
