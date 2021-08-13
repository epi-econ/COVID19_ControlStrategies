#!/bin/bash

#SBATCH --mail-user=abento@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=0-12:0:00
#SBATCH --mem=58gb
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=VSL_5e06

module load r

for VSL in 5e06
do
    for gridsize in 25
    do
	label="VSL_${VSL}"
#	nohup Rscript epicon_script_cheby.R 96 $gridsize 0 "nosocial" "decentralized" 5.166426 7.513051 3.548544 2.6 5.1 58000 0.1 1 "eqm_$label" "R0" 1 0.8554632 1 1 $VSL "VSL" > logs/eqm_$label.log 2>&1

	nohup Rscript epicon_script_cheby.R 48 $gridsize 1 "nosocial" "planner" 5.166426 7.513051 3.548544 2.6 5.1 58000 0.1 1 "planner_$label" "R0" 1.01 0.8554632 1 1 $VSL "VSL" > logs/planner_$label.log 2>&1
    done
done 
