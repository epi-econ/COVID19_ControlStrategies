#!/bin/bash

#SBATCH --mail-user=abento@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=0-12:0:00
#SBATCH --mem=58gb
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=RA_const_VSL_0_1

module load r

for RA in 0.1
do
    for gridsize in 20
    do
	label="const_VSL_RA_${RA}"
	nohup Rscript epicon_script_cheby.R 48 $gridsize 0 "nosocial" "decentralized" 5.166426 7.513051 3.548544 2.6 5.1 58000 $RA 0 "eqm_$label" "R0" 1.01 0.8554632 1 1 "1e07" "VSL" > logs/eqm_$label.log 2>&1

	nohup Rscript epicon_script_cheby.R 48 $gridsize 0 "nosocial" "planner" 5.166426 7.513051 3.548544 2.6 5.1 58000 $RA 0 "planner_$label" "R0" 4 0.8554632 1 1 "1e07" "VSL" > logs/planner_$label.log 2>&1
    done
done 
