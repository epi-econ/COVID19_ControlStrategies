#!/bin/bash 
# shell script to run scenarios, passing scriptargs to epicon_script.R

##### order of arguments to shell_epicon_script:
# ncores = as.numeric(scriptargs[1])
# S,I,R_gridlength = as.numeric(scriptargs[2])
# eta = as.numeric(scriptargs[3])
# hot_start = as.numeric(scriptargs[4])
# scenario_label = as.numeric(scriptargs[5])
# policy = as.numeric(scriptargs[6])
# social = as.character(scriptargs[7])
# high_precision = as.numeric(scriptargs[8])
#####

###########################

##### Check using next-generation method contacts

nohup Rscript epicon_script.R 80 24 0 "nosocial" "planner" 5.166426 7.513051 3.548544 2.6 5.1 58000 0.1 1 "planner_cheby_test_24" "R0" 3 0.8554632 1 1 1e07 "VSL" "cheby"  > planner.log 2>&1 

nohup Rscript epicon_script.R 80 24 0 "nosocial" "decentralized" 5.166426 7.513051 3.548544 2.6 5.1 58000 0.1 1 "eqm_cheby_test_24" "R0" 3 0.8554632 1 1 1e07 "VSL" "cheby"  > eqm.log 2>&1 
